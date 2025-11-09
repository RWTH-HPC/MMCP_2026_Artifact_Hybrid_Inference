// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredbndrycnd.h"

using namespace std;

template <MInt nDim>
StructuredBndryCnd<nDim>::StructuredBndryCnd(FvStructuredSolver<nDim>* solver, StructuredGrid<nDim>* grid)
  : m_StructuredComm(solver->m_StructuredComm),
    m_nCells(solver->m_nCells),
    m_nPoints(solver->m_nPoints),
    m_cells(solver->m_cells),
    CV(solver->CV),
    PV(solver->PV),
    FQ(solver->FQ),
    m_solverId(solver->m_solverId),
    m_noGhostLayers(solver->m_noGhostLayers),
    m_physicalBCMap(solver->m_windowInfo->physicalBCMap),
    m_auxDataMap(solver->m_windowInfo->physicalAuxDataMap),
    m_globalStructuredBndryMaps(solver->m_windowInfo->globalStructuredBndryCndMaps),
    m_noSpongeDomainInfos(solver->m_noSpongeDomainInfos),
    m_spongeBcWindowInfo(solver->m_spongeBcWindowInfo),
    m_spongeLayerType(solver->m_spongeLayerType),
    m_spongeLayerThickness(solver->m_spongeLayerThickness),
    m_sigmaSponge(solver->m_sigmaSponge),
    m_betaSponge(solver->m_betaSponge),
    m_targetDensityFactor(solver->m_targetDensityFactor),
    m_noCells(grid->m_noCells),
    m_sutherlandPlusOne(solver->m_sutherlandPlusOne),
    m_sutherlandConstant(solver->m_sutherlandConstant) {
  TRACE();

  m_solver = solver;
  m_grid = grid;

  // the surfaces for the channel flow calculation
  m_channelSurfaceIn = F0;
  m_channelSurfaceOut = F0;
}


template <MInt nDim>
void StructuredBndryCnd<nDim>::applyDirichletNeumannBC() {
  // TODO_SS labels:FV,totest check if it is eventually necessary to first loop over the 6000er BCs
  for(MInt bcId = 0; bcId < m_noBndryCndIds; bcId++) {
    (this->*bndryCndHandler[bcId])(bcId);
  }
}

template <MInt nDim>
void StructuredBndryCnd<nDim>::applyNonReflectingBC() {
  TRACE();
  // not implemented yet
  // for( MInt bcId = 0;  bcId < m_noBndryCndIds;  bcId++ )
  //  (this->*nonReflectingBoundaryCondition[bcId]) (bcId);
}


template <MInt nDim>
void StructuredBndryCnd<nDim>::assignBndryCnds() {
  TRACE();

  m_noBndryCndIds = m_physicalBCMap.size();

  bndryCndHandler = nullptr;
  bndryCndHandler = new BndryCndHandler[m_physicalBCMap.size()];
  initBndryCndHandler = new BndryCndHandler[m_physicalBCMap.size()];

  // relation between surface index map and the bcId
  mAlloc(m_channelSurfaceIndexMap, std::max((MInt)m_physicalBCMap.size(), 1), "m_channelSurfaceIndexMap", -1, AT_);
  // plenumSurface relation between surface/map index and the bcId
  mAlloc(m_plenumSurfaceIndexMap, std::max((MInt)m_physicalBCMap.size(), 1), "m_plenumSurfaceIndexMap", -1, AT_);
  MInt counter = 0, counterPlenum = 0;
  // assign the function pointers
  for(MUint bcId = 0; bcId < m_physicalBCMap.size(); bcId++) {
    switch(m_physicalBCMap[bcId]->BC) {
      case 0: {
        // euler wall
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc0;         // empty boundary condition
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc0; // empty boundary condition
        break;
      }
      case 1000:
      case 1004: {
        // adiabatic wall
        if(m_solver->m_movingGrid) {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc1004;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc1004;
        } else {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc1000;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc1000;
        }
        break;
      }
      case 1003:
      case 1006: {
        // isothermal wall
        if(m_solver->m_movingGrid) {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc1006;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc1006;
        } else {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc1003;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc1003;
        }
        break;
      }
      case 1001: {
        // euler wall
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc1001;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc1001;
        break;
      }
      case 1007: {
        // oscillating wall
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc1007;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc1007;
        break;
      }
      case 2001: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2001;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2001;
        break;
      }
      case 2002:
      case 2010: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2002;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2002;
        break;
      }
      case 2003: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2003;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2003;
        break;
      }
      case 2004:
      case 2024: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2004;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2004;
        break;
      }
      case 2005: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2005;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2005;
        break;
      }
      case 2006: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2006;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2006;
        break;
      }
      case 2007: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2007;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2007;
        break;
      }
      case 2009: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2009;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2009;
        break;
      }
      case 2012: { // characteristic inflow
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2012;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2012;
        break;
      }
      case 2013: { // characteristic outflow
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2013;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2013;
        break;
      }
      case 2014: { // characteristic outflow
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2014;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2014;
        break;
      }
      case 2020: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2020;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2020;
        break;
      }
      case 2021: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2021;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2021;
        break;
      }
      case 2097: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2097;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2097;
        m_plenumSurfaceIndexMap[bcId] = counterPlenum;
        counterPlenum++;
        break;
      }
      case 2099: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2099;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2099;
        break;
      }
      case 2221: { // junoh //zonal with STG
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2221;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2221;
        break;
      }
      case 2222: { // junoh //zonal without STG
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2222;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2222;
        break;
      }
      case 2199: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2199;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2199;
        break;
      }
      case 2401:
      case 2402: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2402;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2402;
        m_channelSurfaceIndexMap[bcId] = counter;
        counter++;
        break;
      }
      case 2500: { // Rescaling
        if(m_solver->m_rans) {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc2510;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2510;
        } else {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc2500;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2500;
        }
        break;
      }
      case 2501: { // Rescaling
        if(m_solver->m_rans) {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc2511;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2501;
        } else {
          bndryCndHandler[bcId] = &StructuredBndryCnd::bc2501;
          initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2501;
        }
        break;
      }
      case 2600: { // Prescribing profile
        m_solver->m_bc2600 = true;
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2600;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2600;
        m_solver->m_bc2600 = true;
        break;
      }
      case 2601: { // Prescribing profile
        m_solver->m_bc2601 = true;
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2601;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2601;
        m_solver->m_bc2601 = true;
        break;
      }
      case 2700: { // mode inflow
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2700;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2700;
        break;
      }
      case 2730: { // fsc outflow
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2730;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2730;
        break;
      }
      case 2888: { // fsc inflow
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2888;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2888;
        break;
      }
      case 2999: { // Blasius inflow
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2999;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2999;
        break;
      }
      case 2900: { // Jet Freund Inlet
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc2900;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc2900;
        break;
      }
      case 3000: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc3000;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc3000;
        break;
      }
      case 3001: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc3001;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc3001;
        break;
      }
      case 6000: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc6000;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc6000;
        break;
      }
      // case 6002: {
      //   // Fluid-porous interface
      //   bndryCndHandler[bcId] = &StructuredBndryCnd::bc6002;
      //   initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc6002;
      //   break;
      // }
      case 7909: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc7909;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc7909;
        break;
      }
        // empty BC!!!!!!!!!!!!!!!!
      case 4401:
      case 4402:
      case 4403:
      case 4404:
      case 4405:
      case 4406: {
        bndryCndHandler[bcId] = &StructuredBndryCnd::bc9999;
        initBndryCndHandler[bcId] = &StructuredBndryCnd::initBc9999;
        break;
      }

      default: {
        cout << "boundary condtition is missing" << m_physicalBCMap[bcId]->BC << endl;
        mTerm(1, AT_, "Boundary Condition is not implemented");
        break;
      }
    }
  }
}

template <MInt nDim>
void StructuredBndryCnd<nDim>::correctBndryCndIndices() {
  // in correcting cell Information
  for(MInt bcId = 0; bcId < m_noBndryCndIds; bcId++) {
    (this->*initBndryCndHandler[bcId])(bcId);
  }
}

template <MInt nDim>
StructuredBndryCnd<nDim>::~StructuredBndryCnd() {
  delete[] bndryCndHandler;
  delete[] initBndryCndHandler;
}

template <MInt nDim>
void StructuredBndryCnd<nDim>::saveAuxData() {}

// Explicit instantiations for 2D and 3D
template class StructuredBndryCnd<2>;
template class StructuredBndryCnd<3>;
