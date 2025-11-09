// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDBNDRYCND2D
#define STRUCTUREDBNDRYCND2D

#include "fvstructuredbndrycnd.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver2d.h"
#include "fvstructuredwindowmapping.h"

class FvStructuredSolver2D;

/** \brief Class for the 2D stuctured boundary conditions
 *
 */
template <MBool isRans>
class StructuredBndryCnd2D : public StructuredBndryCnd<2> {
 public:
  friend class FvStructuredSolver2D;

  using Timers = maia::structured::Timers_;
  StructuredBndryCnd2D(FvStructuredSolver<2>* solver, StructuredGrid<2>* grid);


  ~StructuredBndryCnd2D();

  static constexpr const MInt nDim = 2;

  void correctBndryCndIndices();
  void correctWallDistanceAtBoundary(MInt);

  void bc1000(MInt); // wall no slip
  template <RansMethod ransMethod>
  void bc1000_(MInt);
  void bc1001(MInt);                  // euler wall
  void bc1003(MInt);                  // isothermal wall
  void bc1004(MInt);                  // moving adiabatic wall
  void bc2001(MInt);                  // subsonic inflow
  void bc2003(MInt);                  // subsonic outflow
  void bc2004(MInt);                  // subsonic outflow
  void bc2002(MInt);                  // supersonic inflow
  void bc2005(MInt);                  // supersonic outflow
  void bc2006(MInt);                  // subsonic zero velocity in/outflow
  void bc2007(MInt);                  // subsonic outflow
  void bc2021(MInt);                  // shear flow inflow
  void bc2199(MInt);                  // subsonic inflow compressible bernoulli (tfs2099)
  void bc2402(MInt);                  // channel flow
  void bc2510(MInt);                  // rescaling inlet
  void bc2511(MInt);                  // rescaling recycling
  void bc2600(MInt);                  // prescribing
  void bc2999(MInt);                  // blasius inflow
  void bc3000(MInt);                  // symmetrie
  virtual void bc6002(MInt) override; // Fluid-porous interface

  void initBc1000(MInt);
  void initBc1001(MInt){};
  void initBc1003(MInt);
  void initBc1004(MInt);
  void initBc2001(MInt);
  void initBc2002(MInt);
  void initBc2003(MInt){};
  void initBc2004(MInt);
  void initBc2005(MInt);
  void initBc2006(MInt);
  void initBc2007(MInt);
  void initBc2021(MInt);
  void initBc2199(MInt){};
  void initBc2402(MInt); // channel flow
  void initBc2500(MInt){};
  void initBc2501(MInt){};
  void initBc2510(MInt);
  void initBc2511(MInt){};
  void initBc2600(MInt);
  void initBc3000(MInt);
  void initBc6002(MInt) override{}; // Fluid-porous interface

  virtual void computeFrictionPressureCoef(MBool computePower) override {
    computeFrictionPressureCoef_<true, true, true>(false, computePower);
  }
  virtual void computeFrictionCoef() override { computeFrictionPressureCoef_<false, true, false>(false, false); };
  template <MBool calcCp, MBool calcCf, MBool calcIntegrals>
  void computeFrictionPressureCoef_(const MBool auxDataWindow = false, const MBool computePower = false);
  template <MBool calcCp, MBool calcCf, MBool interface>
  void calc_cp_cf(const MInt, const MInt, const MInt, const MInt, MFloat (&)[calcCp + nDim * calcCf]);
  virtual void distributeWallAndFPProperties() override;
  void distributeMapProperties(const std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>&,
                               const std::vector<MInt>&,
                               const std::vector<MInt>&,
                               const std::map<MInt, std::tuple<MInt, MInt, MFloat>>& cellId2recvCell,
                               const std::vector<MInt>&,
                               MFloat* const);

  void readAndDistributeSpongeCoordinates();
  void updateSpongeLayer();
  void computeWallDistances();
  virtual void computeLocalWallDistances() override;
  void computeDistance2Map(const std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>&,
                           MFloat* const,
                           std::vector<std::pair<MInt, MInt>>&,
                           std::vector<MInt>&,
                           std::vector<MFloat>&);

  // Placeholder cellId can be used for more sophisticated comparisions
  template <typename T>
  struct comp {
    MBool operator()(const MInt /*cellId*/, const T& a, const T& b) { return a < b; }
  };
  template <typename T = comp<MFloat>>
  void getCloserMap(const MFloat* const,
                    std::vector<std::pair<MInt, MInt>>&,
                    const MFloat* const,
                    std::vector<std::pair<MInt, MInt>>&,
                    MFloat* const,
                    T comparator = {});
  void setUpNearMapComm(const std::vector<std::pair<MInt, MInt>>&,
                        const std::vector<MInt>&,
                        const std::vector<MFloat>&,
                        std::map<MInt, std::tuple<MInt, MInt, MFloat>>&,
                        std::vector<MInt>&,
                        std::vector<MInt>&);
  virtual void computeLocalExtendedDistancesAndSetComm() override;
  void modifyFPDistance(const std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>&,
                        MFloat* const,
                        std::vector<std::pair<MInt, MInt>>&,
                        const std::vector<MFloat>&);
  MFloat shortestDistanceToLineElement(const MFloat (&)[nDim], const MFloat (&)[nDim], const MFloat (&)[nDim], MFloat&,
                                       MFloat&);
  inline MInt cellIndex(MInt i, MInt j);
  inline MInt getPointIdFromCell(MInt i, MInt j);
  inline MInt getPointIdFromPoint(MInt origin, MInt incI, MInt incJ);
  inline MFloat pressure(MInt);
  inline MFloat temperature(MInt);

 protected:
  FvStructuredSolver2D* m_solver;
  MInt m_startCommPeriodic;
  MInt m_endCommPeriodic;
  MInt m_startCommChannel;
  MInt m_endCommChannel;
  MInt m_channelInflowRank;

  // periodic
  MInt m_noPeriodicConnections;

  // 2500 rescaling bc
  MFloat m_rescalingBLT;
  MFloat m_isothermalWallTemperature;
  MFloat m_bc2021Gradient;
};


#endif
