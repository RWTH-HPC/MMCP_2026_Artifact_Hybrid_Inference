// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBINTERFACE_H
#define LBINTERFACE_H
#include <iostream>
#include <vector>

#include "enums.h"
#include "variables.h"

template <class T>
class Collector;
template <MInt nDim>
class LbSolver;
class LbInterfaceCell;
class LbParentCell;
template <MInt nDim>
class CartesianGrid;

/** Base class for the treatment of refinement for the LB module
 *
 * Implements the data structure and operations for the treatment
 * of the interface between cells of different level.
 *
 */
template <MInt nDim>
class LbInterface {
 public:
  template <MInt nDim_>
  friend class LbSolver;

 public:
  LbInterface(LbSolver<nDim>* solver);
  virtual ~LbInterface();

  // Main adaptation functions
  virtual void refineCell(const MInt parentId, const MInt* childIds) = 0;
  virtual void removeChildren(const MInt parentId) = 0;

  // Debug Functions:
  virtual void printInterfaceCells();
  /// Sets the interface cells to defined values (to be watched e.g. in DX)
  void colorInterface();

 protected:
  MInt m_interfaceId{};

  MInt m_noDistributions;

  std::vector<Collector<LbInterfaceCell>*> m_interfaceChildren;
  std::vector<Collector<LbParentCell>*> m_interfaceParents;

  // Use an external pressure force
  MBool m_cellDependentForcing;
  MBool m_externalForcing;
  MBool m_isEELiquid;

  MInt m_methodId;

  MFloat* m_Fext{};
  MFloat* m_Fg{};

  MPrimitiveVariables<nDim>* PV;


  // Flow variables for prolongation and restriction
  MFloat m_nu;

  MString m_interfaceMethod;
  MString m_adaptationInitMethod;
  MInt m_noCoefficients{};

  MInt m_isThermal;
  MInt m_innerEnergy;

  LbSolver<nDim>* m_solver;
};

#endif
