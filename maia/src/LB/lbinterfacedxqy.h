// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBINTERFACEDXQY_H
#define LBINTERFACEDXQY_H

#include "lbinterface.h"
#include "lblatticedescriptor.h"
#include "lbsolverdxqy.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy;
/**
 * \brief Interface class holding all relevant data and methods for treating
 * prolongation, restriction and initialization of new cells for the DxQy
 * lattice model
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbInterfaceDxQy : public LbInterface<nDim> {
 public:
  template <MInt nDim_, MInt nDist_, class SysEqn_>
  friend class LbSolverDxQy;

 public:
  // Add fields used from template base class to avoid calling with 'this->'
  using LbInterface<nDim>::m_interfaceChildren;
  using LbInterface<nDim>::m_interfaceParents;
  using LbInterface<nDim>::PV;
  using LbInterface<nDim>::m_cellDependentForcing;
  using LbInterface<nDim>::m_externalForcing;
  using LbInterface<nDim>::m_Fext;
  using LbInterface<nDim>::m_isEELiquid;
  using LbInterface<nDim>::m_Fg;
  using LbInterface<nDim>::m_isThermal;
  using LbInterface<nDim>::m_innerEnergy;
  using LbInterface<nDim>::m_interfaceMethod;
  using LbInterface<nDim>::m_adaptationInitMethod;

  using Ld = LbLatticeDescriptor<nDim, nDist>;

  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;

  LbInterfaceDxQy(LbSolver<nDim>* solver);
  virtual ~LbInterfaceDxQy();

  // Main interface functions
  void prolongation();
  void restriction();

  virtual void refineCell(const MInt parentId, const MInt* childIds) override final;
  virtual void removeChildren(const MInt parentId) override final;

 protected:
  // Interface functions
  virtual void prolongation0(){};
  virtual void restriction0(){};

  virtual void prolongation10();
  virtual void restriction10();

  template <MBool compressible = false>
  void prolongationDupuis_();
  template <MBool compressible = false>
  void restrictionDupuis_();

  virtual void prolongationDupuis();
  virtual void restrictionDupuis();

  void prolongationDupuisCompressible();
  void restrictionDupuisCompressible();

  virtual void prolongationRohde();
  virtual void restrictionRohde();

  virtual void prolongationThermalDupuis();
  virtual void restrictionThermalDupuis();

  virtual void prolongationThermalRohde();
  virtual void restrictionThermalRohde();

  // Functions for adaptive mesh refinement
  virtual void refineCellDupuis(const MInt parentId, const MInt* childIds);
  virtual void refineCellCopyPaste(const MInt parentId, const MInt* childIds);
  virtual void removeChildsDupuisFilippova(const MInt parentId);
  virtual void removeChildsCopyPaste(const MInt parentId);

  inline void getCellForcing(std::array<MFloat, nDist>& F, const MInt cellId = -1);

 private:
  // sets the prolongation and restriction routines
  void setInterfaceFunctions();

  typedef void (LbInterfaceDxQy::*InterfaceFunction)();
  InterfaceFunction fProlongation;
  InterfaceFunction fRestriction;

  // Sets the refine and coarsen functions
  void setAdaptationFunctions();

  typedef void (LbInterfaceDxQy::*RefineCellFunction)(const MInt parentId, const MInt* childIds);
  RefineCellFunction fRefineCell;

  typedef void (LbInterfaceDxQy::*RemoveChildrenFunction)(const MInt parentId);
  RemoveChildrenFunction fRemoveChildren;

  // Data members used to replace local static variables
  // TODO: Remove as soon as prolong/restrict10 is cleaned up or removed..

  // prolongation10
  MFloat m_static_prolongation10_tmp{};
  MFloat m_static_prolongation10_tmp2{};
  MFloat m_static_prolongation10_b[2 * nDim]{};
  MFloat m_static_prolongation10_c[nDim * nDim]{};
  MFloat m_static_prolongation10_trace{};
  MInt m_static_prolongation10_tmpDistId{};
  // restriction10
  MFloat m_static_restriction10_tmp{};
  MFloat m_static_restriction10_tmp2{};
  MFloat m_static_restriction10_b[2 * nDim]{};
  MFloat m_static_restriction10_c[nDim * nDim]{};
  MFloat m_static_restriction10_trace{};
  MInt m_static_restriction10_tmpDistId{};
};

/**
 * \brief General accessor for potentially cell dependent forcing
 *
 * \param[out] F      Forcing term
 * \param[in]  cellId Cell id (-1 if global forcing is needed)
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbInterfaceDxQy<nDim, nDist, SysEqn>::getCellForcing(std::array<MFloat, nDist>& F,
                                                                 const MInt cellId /*= -1*/) {
  if(m_externalForcing) {
    std::copy(m_Fext, m_Fext + nDist, F.begin());
  } else {
#ifdef WAR_NVHPC_PSTL
    for(MInt d = 0; d < nDist; d++) {
      F[d] = 0.0;
    }
#else
    F.fill(0.0);
#endif
  }
  if(cellId < 0) {
    return;
  }
  if(m_isEELiquid) {
    const MFloat alpha = m_solver->a_alphaGasLim(cellId);
    for(MInt dist = 0; dist < nDist; dist++) {
      F[dist] += m_Fg[dist] * alpha;
    }
  }
}

#endif
