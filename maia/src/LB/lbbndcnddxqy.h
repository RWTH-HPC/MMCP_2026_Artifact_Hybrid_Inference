// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBBNDCNDDXQY_H
#define LBBNDCNDDXQY_H

#include <random>
#include "lbbndcnd.h"
#include "lbconstants.h"
#include "lblatticedescriptor.h"
#include "lbsolverdxqy.h"

using namespace lbconstants;

template <MInt nDim, MInt nDist, class SysEqn>
class LbBndCndDxQy : public LbBndCnd<nDim> {
  const SysEqn m_sysEqn{};

 public:
  // Add fields used from template base class to avoid calling with 'this->'
  using LbBndCnd<nDim>::m_solverId;
  using LbBndCnd<nDim>::m_initialNoCells;
  using LbBndCnd<nDim>::m_rho1;
  using LbBndCnd<nDim>::m_rho2;
  using LbBndCnd<nDim>::m_lbControlInflow;
  using LbBndCnd<nDim>::m_bndCells;
  using LbBndCnd<nDim>::m_bndCndIds;
  using LbBndCnd<nDim>::m_rhoLast;
  using LbBndCnd<nDim>::m_lRho;
  using LbBndCnd<nDim>::m_currentTimeStep;
  using LbBndCnd<nDim>::m_bndCndSegIds;
  using LbBndCnd<nDim>::m_mapBndCndIdSegId;
  using LbBndCnd<nDim>::m_bndCndOffsets;
  using LbBndCnd<nDim>::m_interpolationDistMethod;
  using LbBndCnd<nDim>::m_outputWallDistanceField;
  using LbBndCnd<nDim>::m_gridCutTest;
  using LbBndCnd<nDim>::m_noDistributions;
  using LbBndCnd<nDim>::m_noPeriodicSegments;
  using LbBndCnd<nDim>::m_periodicSegmentsIds;
  using LbBndCnd<nDim>::m_noInOutSegments;
  using LbBndCnd<nDim>::m_inOutSegmentsIds;
  using LbBndCnd<nDim>::m_mapSegIdsInOutCnd;
  using LbBndCnd<nDim>::m_Ma;
  using LbBndCnd<nDim>::m_initialVelocityVecs;
  using LbBndCnd<nDim>::m_blasius_delta;
  using LbBndCnd<nDim>::m_blasius;
  using LbBndCnd<nDim>::m_nu;
  using LbBndCnd<nDim>::PV;
  using LbBndCnd<nDim>::m_methodId;
  using LbBndCnd<nDim>::m_pulsatileFrequency;
  using LbBndCnd<nDim>::m_exDirs;
  using LbBndCnd<nDim>::m_exWeights;
  using LbBndCnd<nDim>::m_omega;
  using LbBndCnd<nDim>::m_densityFluctuations;
  using LbBndCnd<nDim>::m_domainLength;
  using LbBndCnd<nDim>::m_BCComm;
  using LbBndCnd<nDim>::m_BCWallMBComm;
  using LbBndCnd<nDim>::m_totalNoBcCells;
  using LbBndCnd<nDim>::m_referenceLength;
  using LbBndCnd<nDim>::m_Re;
  using LbBndCnd<nDim>::m_deltaRho;
  using LbBndCnd<nDim>::m_maxDeltaRho;
  using LbBndCnd<nDim>::m_ReLast;
  using LbBndCnd<nDim>::m_BCneighbors;
  using LbBndCnd<nDim>::m_noBCNeighbors;
  using LbBndCnd<nDim>::m_BCOutputFileName;
  using LbBndCnd<nDim>::m_BCWallMBNeighbors;
  using LbBndCnd<nDim>::m_allDomainsHaveBC;
  using LbBndCnd<nDim>::m_BCResidualStream;
  using LbBndCnd<nDim>::m_forceFile;
  using LbBndCnd<nDim>::m_localReCutInterval;
  using LbBndCnd<nDim>::m_localReCutReportInterval;
  using LbBndCnd<nDim>::m_hasLocalReCut;
  using LbBndCnd<nDim>::m_localReCutCells;
  using LbBndCnd<nDim>::m_localReCutDiameter;
  using LbBndCnd<nDim>::m_localReCutRe;
  using LbBndCnd<nDim>::m_localReCutAdpPerc;
  using LbBndCnd<nDim>::m_firstBCinComm;
  using LbBndCnd<nDim>::m_lbNoMovingWalls;
  using LbBndCnd<nDim>::m_segIdMovingWalls;
  using LbBndCnd<nDim>::m_lbWallVelocity;
  using LbBndCnd<nDim>::m_lbNoHeatedWalls;
  using LbBndCnd<nDim>::m_segIdHeatedWalls;
  using LbBndCnd<nDim>::m_lbWallTemperature;
  using LbBndCnd<nDim>::m_noReactivatedCells;
  using LbBndCnd<nDim>::m_boundaryCellsMb;
  using LbBndCnd<nDim>::m_boundaryCellMappingMb;
  using LbBndCnd<nDim>::m_calcWallForces;
  using LbBndCnd<nDim>::m_calcWallForcesInterval;
  using LbBndCnd<nDim>::m_bndNormals;
  using LbBndCnd<nDim>::m_latentHeat;
  using LbBndCnd<nDim>::m_calcSublayerDist;

  // Introduce cell type from parent class
  using Cell = typename LbBndCnd<nDim>::Cell;
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  template <MInt nDim_, MInt nDist_, class SysEqn_>
  friend class LbSolverDxQy;

 private:
  MInt bounceBackSchemeMb;
  void (LbBndCndDxQy<nDim, nDist, SysEqn>::*bounceBackFunctionMb)(const MInt cellIndex, const MInt set);

#ifdef WAR_NVHPC_PSTL
  MFloat** m_distances{};
  MInt m_bounceBackFunctionMbCase = 0;
#endif
  MFloat m_zeroInflowVelocity;
  MInt m_refillMethodOrder;
  std::vector<MInt> noMissDistBnd;
  std::vector<std::vector<MFloat>> noMissDistBndWeighted;
  std::vector<MInt> ibbBndIds;

  // mucosa model
  struct {
    MFloat T_o = F0;
    MFloat C_o = F0;
    MFloat mucosaThickness = F0;
    MFloat condRatio = F0;
    MFloat diffRatio = F0;
    MFloat refDx = F0;
    MFloat refT = F0;
    MFloat refC = F0;
    MFloat refDiff = F0;
    MFloat refCondF = F0;
    MInt plugFlow = 0;
  } m_mucosaModel;

  MFloat** m_oldWallTemp = nullptr;
  MInt** m_distIntersectionElementId = nullptr;
  MFloat** m_mucousDist = nullptr;
  MFloat** m_fluidDist = nullptr;

  void calculateSublayerDistances(MInt index);
  inline void calculateWallInterface(MInt cellId, MFloat* wallConc, MFloat* wallTemp);

 public:
  LbBndCndDxQy(LbSolver<nDim>* solver);
  virtual ~LbBndCndDxQy();

  constexpr SysEqn sysEqn() const { return m_sysEqn; }

 protected:
  void calculateWallDistances();
  void calculateWallDistances2D();
  void addWallDistanceFieldToOutputFile(ParallelIo& parallelIo, const MBool writeHeader,
                                        const MBool writeData) override;
  void initVortices(MInt /*index*/){};
  void writeBoundaryVTK(MFloat** vertices, MInt num, MInt segmentId);
  void parabolicInflow(MInt index);
  MFloat** allOwnersGetBoundaryVertices(MInt segmentId, MInt* num);

  void bc0(MInt index);

  void bc66666(MInt set);
  void bc66668(MInt set);

  void bc20000(MInt index);
  void bc20001(MInt index);
  void bc20002(MInt index);
  void bc20003(MInt index);
  void bc20004(MInt index);
  void bc20005(MInt index);
  void bc20020(MInt index);
  void bc20022(MInt index);
  void bc20023(MInt index);
  void bc20024(MInt index);
  void bc20025(MInt index);
  void bc20026(MInt index);
  void bc20027(MInt index);
  void bc20030(MInt index);
  template <MInt direction>
  void slidingWall(MInt index);
  void bc20050(MInt index) { slidingWall<0>(index); };
  void bc20051(MInt index) { slidingWall<1>(index); };
  void bc20052(MInt index) { slidingWall<2>(index); };
  void bc20053(MInt index) { slidingWall<3>(index); };
  void bc20054(MInt index) { slidingWall<4>(index); };
  void bc20055(MInt index) { slidingWall<5>(index); };

  void bc20220(MInt index);

  void bc20226(MInt index) { heatFluxBc(index, 0); };
  void bc20227(MInt index) { heatFluxBc(index, 1); };
  void bc20228(MInt index) { heatFluxBc(index, 2); };
  void heatFluxBc(MInt index, MInt bcMode);
  void bcIBBNeumannInit(MInt index);

  void bc20230(MInt index);

  void bc20501(MInt index);
  void bc20501_init(MInt index);

  void bc10000(MInt index);
  void bc10001(MInt index);
  void bc10002(MInt index);
  void bc10003(MInt index);
  void bc10004(MInt index);
  void bc10010(MInt index);
  void bc10020(MInt index);
  void bc10022(MInt index);
  void bc40000(MInt index);
  void bc40020(MInt index);
  void dnt(MInt index, MInt direction);
  void bc40030(MInt index);
  void bc40130(MInt index);
  void bc10050(MInt index);
  void bc10060(MInt index);
  void bc10061(MInt index);

  void bc10090(MInt index);

  void bc10070(MInt index);
  void bc10080(MInt index);
  void bc10111(MInt index);
  void bc40070(MInt index);
  void bc40071(MInt index);
  void bc40072(MInt index);
  void recalcRho(MInt index);
  void bc40073(MInt index);
  void bc40080(MInt index);
  void bc40081(MInt index);
  void bc40082(MInt index);
  void bc40072_40082_init(MInt index);

  void bc40100(MInt index);
  void bc40110(MInt index);

  void bc40120(MInt index);

  void bc30000(MInt index);

  template <MUint type>
  inline void calcCharValues(const MInt index, const MInt bndCellId, MFloat& rho_b, MFloat* u_b);
  inline void calcCharValuesOnFace(const MInt index, const MInt direction, const MInt bndCellId, MFloat& rho_b,
                                   MFloat* u_b);
  void charVelocity(MInt index, MInt direction);
  void bc10040(MInt index) { charVelocity(index, 0); };
  void bc10041(MInt index) { charVelocity(index, 1); };
  void bc10042(MInt index) { charVelocity(index, 2); };
  void bc10043(MInt index) { charVelocity(index, 3); };
  void bc10044(MInt index) { charVelocity(index, 4); };
  void bc10045(MInt index) { charVelocity(index, 5); };
  void bc10046(MInt index);

  void charPressure(MInt index, MInt direction);
  void bc40040(MInt index) { charPressure(index, 0); };
  void bc40041(MInt index) { charPressure(index, 1); };
  void bc40042(MInt index) { charPressure(index, 2); };
  void bc40043(MInt index) { charPressure(index, 3); };
  void bc40044(MInt index) { charPressure(index, 4); };
  void bc40045(MInt index) { charPressure(index, 5); };
  void bc40046(MInt index);

  void charPressure2(MInt index, MInt direction);
  void bc40050(MInt index) { charPressure2(index, 0); };
  void bc40051(MInt index) { charPressure2(index, 1); };
  void bc40052(MInt index) { charPressure2(index, 2); };
  void bc40053(MInt index) { charPressure2(index, 3); };
  void bc40054(MInt index) { charPressure2(index, 4); };
  void bc40055(MInt index) { charPressure2(index, 5); };

  template <MInt direction>
  void outflow(MInt index);
  void bc30010(MInt index) { outflow<0>(index); };
  void bc30011(MInt index) { outflow<1>(index); };
  void bc30012(MInt index) { outflow<2>(index); };
  void bc30013(MInt index) { outflow<3>(index); };
  void bc30014(MInt index) { outflow<4>(index); };
  void bc30015(MInt index) { outflow<5>(index); };

  template <MBool thermal>
  void slipFlow(MInt index);
  void bc30020(MInt index) { slipFlow<false>(index); };
  void bc30021(MInt /*index*/) { TERMM(1, "Only use 3020 BC, not 3021"); };
  void bc30022(MInt /*index*/) { TERMM(1, "Only use 3020 BC, not 3022"); };
  void bc30023(MInt /*index*/) { TERMM(1, "Only use 3020 BC, not 3023"); };
  void bc30024(MInt /*index*/) { TERMM(1, "Only use 3020 BC, not 3024"); };
  void bc30025(MInt /*index*/) { TERMM(1, "Only use 3020 BC, not 3025"); };

  template <MInt direction>
  void outflowLinear(MInt index);
  void bc30030(MInt index) { outflowLinear<0>(index); };
  void bc30031(MInt index) { outflowLinear<1>(index); };
  void bc30032(MInt index) { outflowLinear<2>(index); };
  void bc30033(MInt index) { outflowLinear<3>(index); };
  void bc30034(MInt index) { outflowLinear<4>(index); };
  void bc30035(MInt index) { outflowLinear<5>(index); };

  void bc30040(MInt index) { slipFlow<true>(index); };
  void bc30041(MInt /*index*/) { TERMM(1, "Only use 3040 BC, not 3041"); };
  void bc30042(MInt /*index*/) { TERMM(1, "Only use 3040 BC, not 3042"); };
  void bc30043(MInt /*index*/) { TERMM(1, "Only use 3040 BC, not 3043"); };
  void bc30044(MInt /*index*/) { TERMM(1, "Only use 3040 BC, not 3044"); };
  void bc30045(MInt /*index*/) { TERMM(1, "Only use 3040 BC, not 3045"); };

  void pab(MInt index);
  void bc40060(MInt index) { pab(index); };

  void writeBCOutput(MInt index);

  void interpolatedBounceBackSingleSpecies(const MInt cellId, const MFloat* const uW);
  void interpolatedBounceBackSingleSpeciesThermal(const MInt cellId, const MFloat* const wTPtr, const MFloat* const uW);
  void interpolatedBounceBackSingleSpeciesThermal(const MInt cellId, const MFloat wT, const MFloat* const uW) {
    std::array<MFloat, nDist> wTPtr{};
    for(MInt d = 0; d < nDist; d++) {
      wTPtr[d] = wT;
    }
    interpolatedBounceBackSingleSpeciesThermal(cellId, wTPtr.data(), uW);
  };
  void interpolatedBounceBackSingleSpeciesTransport(const MInt cellId, const MFloat* const wCPtr,
                                                    const MFloat* const uW);
  void interpolatedBounceBackSingleSpeciesTransport(const MInt cellId, const MFloat wC, const MFloat* const uW) {
    std::array<MFloat, nDist> wCPtr{};
    for(MInt d = 0; d < nDist; d++) {
      wCPtr[d] = wC;
    }
    interpolatedBounceBackSingleSpeciesTransport(cellId, wCPtr.data(), uW);
  };
  void interpolatedBounceBackSingleSpeciesThermalFlux(const MInt cellId, const MFloat qT, MInt bcIndex);
  inline void calculateEqDistsWallSingleSpecies(const MInt pCellId);
  inline void calculateEqDistsWallSingleSpeciesThermal(const MInt pCellId, MFloat wT);
  inline void calculateEqDistsWallSingleSpeciesTransport(const MInt pCellId, MFloat wC);
  void calculateWallForces(MInt index);
  void calculateWallForcesMb(MInt set);
  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;
  inline void extrapolateVariable(const MInt index, const MInt pCellId, const MInt var, MFloat* const p_var);
  inline void extrapolateVelocities(MInt index, const MInt pCellId, MFloat* l_uu);
  inline void incrementForces(const MInt cellId, const MInt j, const MFloat* uW, MFloat* forces);
  inline void incrementForces(const MInt cellId, const MInt mbCellId, const MInt j, const MFloat* uW);
  void getBoundaryVelocityMb(const MInt cellId, MFloat* uW);
  void getBoundaryVelocity(const MInt index, MFloat* uW);
  MFloat getBoundaryTemperature(const MInt index);
  inline MFloat firstMomentSourceTerm(const MFloat* const uW, const MFloat rho, const MInt dist);
  template <MInt mode>
  MFloat zerothMomentSourceTerm(const MInt pCellId, const MInt dist, const MFloat var, const MFloat* uW);

  inline void calculateEqDistsTotalEnergy(const MInt pCellId, const MFloat l_rho, const MFloat squared_velocity,
                                          const MFloat l_uu[nDim], const MFloat l_t);
  inline void calculateEqDistsTransport(const MInt pCellId, const MFloat l_rho, const MFloat squared_velocity,
                                        const MFloat l_uu[nDim], MFloat l_c);
  inline void interpolatedBounceBackMb_Bouzidi_lin(const MInt cellIndex, const MInt set);
  inline void interpolatedBounceBackMb_Bouzidi_qua(const MInt cellIndex, const MInt set);
  inline void interpolatedBounceBackMb_Yu_qua(const MInt cellIndex, const MInt set);
  inline void interpolatedBounceBackMb_Bouzidi_lin_thermal(const MInt cellIndex, const MInt set);
  inline void interpolatedBounceBackMb_Bouzidi_lin_transport(const MInt cellIndex, const MInt set);

  inline void interpolatedAntiBounceBackMb_Bouzidi_qua(const MInt cellIndex, const MInt set);

 public:
  void extrapolateVelocitiesMb();
  void refillEmergedCell(const MInt /*pCellId*/);
  void refillEmergedCellNormalExtrapolationLinear(const MInt /*pCellId*/);
  void refillEmergedCellNormalExtrapolationQuadratic(const MInt /*pCellId*/);
  inline MFloat getDistanceMb(const MInt cellId, const MInt mbCellId, const MInt distId) {
    const MFloat cellLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(cellId));

    MFloat q;
    if(distId < Ld::distFld(0)) {
      q = m_boundaryCellsMb.distance(mbCellId, distId) / cellLength;
    } else if(distId < (Ld::distFld(0) + Ld::distFld(1))) {
      q = m_boundaryCellsMb.distance(mbCellId, distId) / (SQRT2 * cellLength);
    } else {
      q = m_boundaryCellsMb.distance(mbCellId, distId) / (SQRT3 * cellLength);
    }
    return q;
  }
};
#endif
