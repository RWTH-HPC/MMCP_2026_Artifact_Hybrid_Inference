// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBSOLVERDXQY_H
#define LBSOLVERDXQY_H

#include <random>
#include "lbfunctions.h"
#include "lbinterfacedxqy.h"
#include "lblatticedescriptor.h"
#include "lbsolver.h"
#include "lbsrctermcontroller.h"

#include "lbsyseqn.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbBndCndDxQy;

template <MInt nDim, MInt nDist, class SysEqn>
class LbInterfaceDxQy;

/** \brief This class represents all LB models
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy : public LbSolver<nDim> {
 public:
  const SysEqn m_sysEqn{};

  using Base = LbSolver<nDim>;
  using GridProxy = typename LbSolver<nDim>::GridProxy;
  using LbBndCnd = LbBndCndDxQy<nDim, nDist, SysEqn>;

  LbSolverDxQy(MInt id, MInt noDistributions, GridProxy& gridProxy_, Geometry<nDim>& geometry_, const MPI_Comm comm);
  ~LbSolverDxQy();

  using Base::a_alphaGasLim;
  using Base::a_coordinate;
  using Base::a_diffusivity;
  using Base::a_distribution;
  using Base::a_distributionThermal;
  using Base::a_distributionTransport;
  using Base::a_externalForces;
  using Base::a_hasNeighbor;
  using Base::a_isActive;
  using Base::a_isBndryGhostCell;
  using Base::a_isHalo;
  using Base::a_isInterfaceChild;
  using Base::a_kappa;
  using Base::a_level;
  using Base::a_levelSetFunctionMB;
  using Base::a_noCells;
  using Base::a_nu;
  using Base::a_nuT;
  using Base::a_oldDistribution;
  using Base::a_oldDistributionThermal;
  using Base::a_oldDistributionTransport;
  using Base::a_oldVariable;
  using Base::a_oldVariables_ptr;
  using Base::a_variable;
  using Base::a_variables_ptr;
  using Base::c_cellLengthAtLevel;
  using Base::c_isLeafCell;
  using Base::c_level;
  using Base::c_neighborId;
  using Base::c_noChildren;
  using Base::c_parentId;
  using Base::centerOfGravity;
  using Base::domainId;
  using Base::exchangeData;
  using Base::getReLambdaAndUrmsInit;
  using Base::grid;
  using Base::initNu;
  using Base::m_activeCellList;
  using Base::m_activeCellListLvlOffset;
  using Base::m_adaptation;
  using Base::m_arraySize;
  using Base::m_bandWidth;
  using Base::m_calculateDissipation;
  using Base::m_cells;
  using Base::m_controlVelocity;
  using Base::m_Cs;
  using Base::m_currentMaxNoCells;
  using Base::m_deltaX;
  using Base::m_densityFluctuations;
  using Base::m_densityGradient;
  using Base::m_diffusivity;
  using Base::m_domainLength;
  using Base::m_EELiquid;
  using Base::m_extractedCells;
  using Base::m_Fext;
  using Base::m_FftInit;
  using Base::m_finalRe;
  using Base::m_Ga;
  using Base::m_geometry;
  using Base::m_geometryIntersection;
  using Base::m_gridPoints;
  using Base::m_initDensityGradient;
  using Base::m_initFromCoarse;
  using Base::m_initialNoCells;
  using Base::m_initMethod;
  using Base::m_initRe;
  using Base::m_initStartTime;
  using Base::m_initTime;
  using Base::m_innerEnergy;
  using Base::m_isCompressible;
  using Base::m_isEELiquid;
  using Base::m_isRefined;
  using Base::m_isThermal;
  using Base::m_isTransport;
  using Base::m_kappa;
  using Base::m_levelSetId;
  using Base::m_Ma;
  using Base::m_maxNoSets;
  using Base::m_maxResId;
  using Base::m_MijLij;
  using Base::m_MijMij;
  using Base::m_momentumFlux;
  using Base::m_noEmbeddedBodies;
  using Base::m_noPeakModes;
  using Base::m_nu;
  using Base::m_omega;
  using Base::m_particleMomentumCoupling;
  using Base::m_Pe;
  using Base::m_Pr;
  using Base::m_Re;
  using Base::m_referenceLength;
  using Base::m_ReLambda;
  using Base::m_rescoordinates;
  using Base::m_residual;
  using Base::m_residualInterval;
  using Base::m_resRePos;
  using Base::m_resTimestepPos;
  using Base::m_ReTau;
  using Base::m_solutionInterval;
  using Base::m_tanhInit;
  using Base::m_tanhScaleFactor;
  using Base::m_tmpResidual;
  using Base::m_tmpResidualLvl;
  using Base::m_totalEnergy;
#ifdef WAR_NVHPC_PSTL
  using Base::m_distFld;
  using Base::m_distType;
  using Base::m_faculty;
  using Base::m_idFld;
  using Base::m_mFld1;
  using Base::m_mFld2;
  using Base::m_nFld;
  using Base::m_oppositeDist;
  using Base::m_pFld;
  using Base::m_ppdfDir;
  using Base::m_tp;
#endif
  using Base::m_updateAfterPropagation;
  using Base::m_updateMacroscopicLocation;
  using Base::m_UrmsInit;
  using Base::m_velocityControl;
  using Base::m_volumeAccel;
  using Base::m_volumeAccelBase;
  using Base::maxLevel;
  using Base::minLevel;
  using Base::mpiComm;
  using Base::noDomains;
  using Base::noInternalCells;
  using Base::noVariables;
  using Base::outputDir;
  using Base::PV;
  using Base::saveUVWRhoTOnlyPar;
  using Base::solverMethod;
  using Base::swap_variables;

  using Base::computeFFTStatistics;
  using Base::m_resFileName;
  using Base::mRes;

  using Ld = LbLatticeDescriptor<nDim, nDist>;
  static constexpr MInt m_noDistributions = nDist;

  template <MInt nDim_, MInt nDist_, class SysEqn_>
  friend class LbBndCndDxQy;

  template <MInt nDim_, MInt nDist_, class SysEqn_>
  friend class LbInterfaceDxQy;

  template <MInt nDim_>
  friend class LbInterface;

  template <MInt nDim_, MInt nDist_, class SysEqn_>
  friend class CouplingLB;

  constexpr SysEqn sysEqn() const { return m_sysEqn; }

  virtual void initializeLatticeBgk();
  inline void initLatticeBgkFftPipe();
  inline void initLatticeBgkFftChannel();
  inline void initLatticeBgkFftMixing();
  inline void initLatticeBgkFftMixingFilter();
  inline void initLatticeBgkFftIsotropicTurbulence();
  virtual void initLatticeBgkLaminarChannel();
  virtual void initLatticeBgkTurbulentChannel();
  virtual void initLatticeBgkTurbulentMixing();
  virtual void initLatticeBgkTurbulentBoundary();
  virtual void initLatticeBgkTurbulentPipe();
  virtual void initLatticeBgkTurbulentDuct();
  virtual void initLatticeBgkLaminarPipe();
  virtual void initLatticeBgkVortex();
  virtual void initLatticeBgkTGV();
  void initLatticeBgkSpinningVorticies();
  void initLatticeBgkGaussPulse();
  void initLatticeBgkGaussDiffusion();
  void initLatticeBgkGaussAdvection();
  void initLatticeBgkLaminarCylinder();

  // generalized direction init
  virtual void initLatticeBgkLaminarDir(MInt dir);

  virtual void initLatticeBgkLaminar();

  // General init for PPDFs
  void initEqDistFunctions();
  void initNonEqDistFunctions();
  virtual void initThermalEqDistFunctions();
  virtual void initTransportEqDistFunctions();
  virtual void restartInitLb();
  virtual void propagation_step() override;
  virtual void propagation_step_vol() override;
  virtual void propagation_step_thermal() override;
  virtual void propagation_step_thermal_vol() override;
  virtual void propagation_step_transport() override;
  virtual void propagation_step_transport_vol() override;
  virtual void propagation_step_thermaltransport() override;
  virtual void propagation_step_thermaltransport_vol() override;
  template <MInt timeStepOffset = 0>
  void updateVariablesFromOldDist_();
  virtual void updateVariablesFromOldDist() override;
  virtual void updateVariablesFromOldDist_preCollision() override;
  virtual void volumeForces();
  virtual void controlVelocity();
  virtual void averageGlobalVelocity(const MInt dir);

  // Member functions needed for generic call of functions for adaptation
  virtual void removeChildsLb(const MInt parentId);
  virtual void refineCellLb(const MInt parentId, const MInt* childIds);
  virtual void restartBndCnd();

  void sensorVorticity(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                       std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorDivergence(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                        std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorTotalPressure(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                           std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorInterface(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                       std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorMeanStress(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                        std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;

  template <MBool useSmagorinsky>
  void clb_collision_step_base();
  void clb_collision_step() override;
  void clb_smagorinsky_collision_step() override;
  void cumulant_collision_step() override;
  virtual void bgki_collision_step();
  virtual void bgki_thermal_collision_step();
  virtual void bgki_innerEnergy_collision_step();
  virtual void bgki_totalEnergy_collision_step();
  virtual void bgkc_transport_collision_step();
  virtual void bgkc_thermal_transport_collision_step();
  virtual void bgkc_innerenergy_transport_collision_step();
  virtual void bgkc_totalenergy_transport_collision_step();
  void bgki_collision_step_Guo_forcing() override;
  template <MBool compressible = false>
  inline void calculateMomentumFlux(const MInt pCellId);
  template <MBool compressible = false>
  inline void calculateMomentumFlux(const MInt pCellId, MFloat* const momentumFlux);
  template <MBool compressible = false>
  inline void calculateMacroscopicVariables(const MInt cellId, MFloat& rho, MFloat* const u);
  void calculateDissipation();
  template <MBool optimized, MBool useSmagorinsky>
  void mrt_collision_step_base();
  void mrt_collision_step() override;
  void mrt2_collision_step() override;
  void mrt_smagorinsky_collision_step() override;
  void bgki_euler_collision_step() override;
  void bgki_init_collision_step() override;
  void bgkc_collision_step() override;
  void bgki_dynamic_smago_collision_step() override;
  template <MInt mode>
  void bgki_smagorinsky_collision_step_base();
  template <MInt thermalMode>
  void bgki_thermal_collision_step_base();
  template <MInt thermalMode>
  void bgki_thermal_and_transport_collision_step_base();
  void bgki_smagorinsky_collision_step() override;
  void bgki_smagorinsky_collision_step2() override;
  void bgki_smago_wall_collision_step() override;

  template <MBool useSmagorinsky>
  void rbgk_collision_step_base();
  void rbgk_collision_step() override;
  void rbgk_smagorinsky_collision_step() override;
  void rbgk_dynamic_smago_collision_step() override;

  virtual void calculateResidual();
  void initSrcTermController() override;
  void initSrcTerms() override;
  void preCollisionSrcTerm() override;
  void postCollisionSrcTerm() override;
  void postPropagationSrcTerm() override;
  void postCollisionBc() override;
  void postPropagationBc() override;
  virtual void calcNodalLsValues();
  //  virtual void writeToResidualFile(const MChar* text);
  virtual MBool maxResidual();
  void averageSGSTensors(const MInt direction, MInt& count, std::vector<MFloat>& meanTensors);
  void updateMacroscopicVariables();
  void updateViscosity() override;
  void calculateSGSTensors();

  std::mt19937 randNumGen;
  std::uniform_real_distribution<> distrib{0.0, 1.0};

#ifdef WAR_NVHPC_PSTL
  void initArraysForPSTL();
#endif

 protected:
  //! Pointers for the Boundary Conditions, for flow solving

  LbBndCnd* m_bndCnd;
  maia::lb::LbSrcTermController<nDim, nDist, SysEqn> m_srcTermController;
  LbInterfaceDxQy<nDim, nDist, SysEqn>* m_interface;

 protected:
  //! Pointers for the Boundary Conditions, for flow solving

  void (LbSolverDxQy<nDim, nDist, SysEqn>::*m_initMethodPtr)();

 private:
  void prolongation();
  void restriction();

  void initPressureForce() override;
  void initVolumeForces() override;
  template <MBool compressible>
  void initRunCorrection_();
  void initRunCorrection() override;

  MFloat m_smallestCellLength{};

  // Constants for Law of the Wall
  const MFloat C1, C2, C3, C4;

 public:
  inline void setEqDists(const MInt cellId, const MFloat rho, const MFloat* const velocity);
  inline void setEqDists(const MInt cellId, const MFloat rho, const MFloat squaredVelocity,
                         const MFloat* const velocity);
  inline void setEqDistsThermal(const MInt cellId, const MFloat T, const MFloat rho, const MFloat* const velocity);
  inline void setEqDistsThermal(const MInt cellId, const MFloat T, const MFloat rho, const MFloat squaredVelocity,
                                const MFloat* const velocity);
  template <MInt thermalMode>
  inline void setEqDistsThermal_(const MInt cellId, const MFloat T, const MFloat rho, const MFloat* const velocity);
  template <MInt thermalMode>
  inline void setEqDistsThermal_(const MInt cellId, const MFloat T, const MFloat rho, const MFloat squaredVelocity,
                                 const MFloat* const velocity);
  inline void setEqDistsTransport(const MInt cellId, const MFloat C, const MFloat* const velocity);
  inline void setEqDistsTransport(const MInt cellId, const MFloat C, const MFloat squaredVelocity,
                                  const MFloat* const velocity);

  std::array<MFloat, nDist> getEqDists(const MFloat rho, const MFloat* const velocity);
  std::array<MFloat, nDist> getEqDists(const MFloat rho, const MFloat squaredVelocity, const MFloat* const velocity);
  std::array<MFloat, nDist> getEqDistsThermal(const MFloat t, const MFloat rho, const MFloat* const velocity);
  std::array<MFloat, nDist> getEqDistsThermal(const MFloat t, const MFloat rho, const MFloat squaredVelocity,
                                              const MFloat* const velocity);
  template <MInt thermalMode>
  std::array<MFloat, nDist> getEqDistsThermal_(const MFloat t, const MFloat rho, const MFloat* const velocity);
  template <MInt thermalMode>
  std::array<MFloat, nDist> getEqDistsThermal_(const MFloat t, const MFloat rho, const MFloat squaredVelocity,
                                               const MFloat* const velocity);
  std::array<MFloat, nDist> getEqDistsTransport(const MFloat c, const MFloat* const velocity);
  std::array<MFloat, nDist> getEqDistsTransport(const MFloat c, const MFloat squaredVelocity,
                                                const MFloat* const velocity);
};

// TODO: Unify different equations/systems of equations

/**
 * \brief Set BOTH distributions to equilibrium
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \tparam compressible Determines which equilibrium formulation is used
 *
 * \param[in] cellId   Cell id in which the equilibrium is set
 * \param[in] rho      Macroscropic density
 * \param[in] velocity Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDists(const MInt cellId, const MFloat rho,
                                                          const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  sysEqn().calcEqDists(rho, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(), m_distFld.data());
  for(MInt dist = 0; dist < nDist; dist++) {
    a_distribution(cellId, dist) = eqDist[dist];
    a_oldDistribution(cellId, dist) = eqDist[dist];
  }
#else
  sysEqn().calcEqDists(rho, velocity, eqDist.data());

  std::copy_n(eqDist.begin(), nDist, &a_distribution(cellId, 0));
  std::copy_n(eqDist.begin(), nDist, &a_oldDistribution(cellId, 0));
#endif
}

/**
 * \brief Set BOTH distributions to equilibrium
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \tparam compressible Determines which equilibrium formulation is used
 *
 * \param[in] cellId            Cell id in which the equilibrium is set
 * \param[in] rho               Macroscropic density
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDists(const MInt cellId, const MFloat rho,
                                                          const MFloat squaredVelocity, const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  sysEqn().calcEqDists(rho, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                       m_distFld.data());
  for(MInt dist = 0; dist < nDist; dist++) {
    a_distribution(cellId, dist) = eqDist[dist];
    a_oldDistribution(cellId, dist) = eqDist[dist];
  }
#else
  sysEqn().calcEqDists(rho, squaredVelocity, velocity, eqDist.data());

  std::copy_n(eqDist.begin(), nDist, &a_distribution(cellId, 0));
  std::copy_n(eqDist.begin(), nDist, &a_oldDistribution(cellId, 0));
#endif
}

/**
 * \brief Calls function for setting thermal distributions to equilibrium
 *
 * Use this version if the squaredVelocity is not already precomputed
 *
 * \author Julian Vorspohl, Moritz Waldmann
 *
 * \param[in] cellId            Cell id in which the equilibrium is set
 * \param[in] T                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDistsThermal(const MInt cellId, const MFloat T, const MFloat rho,
                                                                 const MFloat* const velocity) {
  if(m_innerEnergy) {
    setEqDistsThermal_<1>(cellId, T, rho, velocity);
  } else if(m_totalEnergy) {
    setEqDistsThermal_<2>(cellId, T, rho, velocity);
  } else {
    setEqDistsThermal_<0>(cellId, T, rho, velocity);
  }
}

/**
 * \brief Set BOTH thermal distributions to equilibrium
 *
 * Use this version if the squaredVelocity is not already precomputed
 *
 * \author Julian Vorspohl, Moritz Waldmann
 *
 * \tparam thermalMode Determines which equilibrium formulation is used
 *
 * \param[in] cellId            Cell id in which the equilibrium is set
 * \param[in] T                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt thermalMode>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDistsThermal_(const MInt cellId, const MFloat T, const MFloat rho,
                                                                  const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  IF_CONSTEXPR(thermalMode == 0) {
    lbfunc::calcEqDistsThermal<nDim, nDist>(T, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                                            m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 1) {
    lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(T, rho, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                                m_tp.data(), m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 2) {
    lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(T, rho, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                                m_tp.data(), m_distFld.data());
  }
  for(MInt dist = 0; dist < nDist; dist++) {
    a_distributionThermal(cellId, dist) = eqDist[dist];
    a_oldDistributionThermal(cellId, dist) = eqDist[dist];
  }
#else
  IF_CONSTEXPR(thermalMode == 0) { lbfunc::calcEqDistsThermal<nDim, nDist>(T, velocity, eqDist.data()); }
  IF_CONSTEXPR(thermalMode == 1) { lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(T, rho, velocity, eqDist.data()); }
  IF_CONSTEXPR(thermalMode == 2) { lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(T, rho, velocity, eqDist.data()); }
  std::copy_n(eqDist.begin(), nDist, &a_distributionThermal(cellId, 0));
  std::copy_n(eqDist.begin(), nDist, &a_oldDistributionThermal(cellId, 0));
#endif
}

/**
 * \brief Calls function for setting thermal distributions to equilibrium
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Julian Vorspohl, Moritz Waldmann
 *
 * \param[in] cellId            Cell id in which the equilibrium is set
 * \param[in] T                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDistsThermal(const MInt cellId, const MFloat T, const MFloat rho,
                                                                 const MFloat squaredVelocity,
                                                                 const MFloat* const velocity) {
  if(m_innerEnergy) {
    setEqDistsThermal_<1>(cellId, T, rho, squaredVelocity, velocity);
  } else if(m_totalEnergy) {
    setEqDistsThermal_<2>(cellId, T, rho, squaredVelocity, velocity);
  } else {
    setEqDistsThermal_<0>(cellId, T, rho, squaredVelocity, velocity);
  }
}

/**
 * \brief Set BOTH thermal distributions to equilibrium
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Julian Vorspohl, Moritz Waldmann
 *
 * \tparam thermalMode Determines which equilibrium formulation is used
 *
 * \param[in] cellId            Cell id in which the equilibrium is set
 * \param[in] T                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt thermalMode>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDistsThermal_(const MInt cellId, const MFloat T, const MFloat rho,
                                                                  const MFloat squaredVelocity,
                                                                  const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  IF_CONSTEXPR(thermalMode == 0) {
    lbfunc::calcEqDistsThermal<nDim, nDist>(T, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                            m_tp.data(), m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 1) {
    lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(T, rho, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(),
                                                m_mFld2.data(), m_tp.data(), m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 2) {
    lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(T, rho, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(),
                                                m_mFld2.data(), m_tp.data(), m_distFld.data());
  }
  for(MInt dist = 0; dist < nDist; dist++) {
    a_distributionThermal(cellId, dist) = eqDist[dist];
    a_oldDistributionThermal(cellId, dist) = eqDist[dist];
  }
#else
  IF_CONSTEXPR(thermalMode == 0) {
    lbfunc::calcEqDistsThermal<nDim, nDist>(T, squaredVelocity, velocity, eqDist.data());
  }
  IF_CONSTEXPR(thermalMode == 1) {
    lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(T, rho, squaredVelocity, velocity, eqDist.data());
  }
  IF_CONSTEXPR(thermalMode == 2) {
    lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(T, rho, squaredVelocity, velocity, eqDist.data());
  }
  std::copy_n(eqDist.begin(), nDist, &a_distributionThermal(cellId, 0));
  std::copy_n(eqDist.begin(), nDist, &a_oldDistributionThermal(cellId, 0));
#endif
}

/**
 * \brief Set BOTH transport distributions to equilibrium
 *
 * This uses the basic transport equilibrium model
 *
 * \param[in] cellId            Cell id in which the equilibrium is set
 * \param[in] C                 Macroscopic concentration
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDistsTransport(const MInt cellId, const MFloat C,
                                                                   const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  lbfunc::calcEqDistsTransport<nDim, nDist>(C, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                                            m_distFld.data());
  for(MInt dist = 0; dist < nDist; dist++) {
    a_distributionTransport(cellId, dist) = eqDist[dist];
    a_oldDistributionTransport(cellId, dist) = eqDist[dist];
  }
#else
  lbfunc::calcEqDistsTransport<nDim, nDist>(C, velocity, eqDist.data());

  std::copy_n(eqDist.begin(), nDist, &a_distributionTransport(cellId, 0));
  std::copy_n(eqDist.begin(), nDist, &a_oldDistributionTransport(cellId, 0));
#endif
}

/**
 * \brief Set BOTH transport distributions to equilibrium
 *
 * This uses the basic transport equilibrium model
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \param[in] cellId            Cell id in which the equilibrium is set
 * \param[in] C                 Macroscopic concentration
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::setEqDistsTransport(const MInt cellId, const MFloat C,
                                                                   const MFloat squaredVelocity,
                                                                   const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  lbfunc::calcEqDistsTransport<nDim, nDist>(C, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                            m_tp.data(), m_distFld.data());
  for(MInt dist = 0; dist < nDist; dist++) {
    a_distributionTransport(cellId, dist) = eqDist[dist];
    a_oldDistributionTransport(cellId, dist) = eqDist[dist];
  }
#else
  lbfunc::calcEqDistsTransport<nDim, nDist>(C, squaredVelocity, velocity, eqDist.data());

  std::copy_n(eqDist.begin(), nDist, &a_distributionTransport(cellId, 0));
  std::copy_n(eqDist.begin(), nDist, &a_oldDistributionTransport(cellId, 0));
#endif
}

/**
 * \brief Calls function to return the equilibrium distribution
 *
 * Use this version if the squaredVelocity is not already precomputed
 *
 * \author Moritz Waldmann
 *
 * \param[in] rho               Macroscopic density
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDists(const MFloat rho,
                                                                        const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist{};
#ifdef WAR_NVHPC_PSTL
  sysEqn().calcEqDists(rho, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(), m_distFld.data());
#else
  sysEqn().calcEqDists(rho, velocity, eqDist.data());
#endif
  return eqDist;
}

/**
 * \brief Calls function to return the equilibrium distribution
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Moritz Waldmann
 *
 * \param[in] rho               Macroscopic density
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDists(const MFloat rho, const MFloat squaredVelocity,
                                                                        const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist{};
#ifdef WAR_NVHPC_PSTL
  sysEqn().calcEqDists(rho, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                       m_distFld.data());
#else
  sysEqn().calcEqDists(rho, squaredVelocity, velocity, eqDist.data());
#endif
  return eqDist;
}

/**
 * \brief Calls function to return the thermal equilibrium distribution
 *
 * Use this version if the squaredVelocity is not already precomputed
 *
 * \author Moritz Waldmann
 *
 * \param[in] t                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDistsThermal(const MFloat t, const MFloat rho,
                                                                               const MFloat* const velocity) {
  if(m_innerEnergy) {
    return getEqDistsThermal_<1>(t, rho, velocity);
  } else if(m_totalEnergy) {
    return getEqDistsThermal_<2>(t, rho, velocity);
  } else {
    return getEqDistsThermal_<0>(t, rho, velocity);
  }
}

/**
 * \brief Return thermal distributions to equilibrium
 *
 * Use this version if the squaredVelocity is not already precomputed
 *
 * \author Moritz Waldmann
 *
 * \tparam thermalMode Determines which equilibrium formulation is used
 *
 * \param[in] t                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt thermalMode>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDistsThermal_(const MFloat t, const MFloat rho,
                                                                                const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  IF_CONSTEXPR(thermalMode == 0) {
    lbfunc::calcEqDistsThermal<nDim, nDist>(t, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                                            m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 1) {
    lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(t, rho, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                                m_tp.data(), m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 2) {
    lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(t, rho, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                                m_tp.data(), m_distFld.data());
  }
#else
  IF_CONSTEXPR(thermalMode == 0) { lbfunc::calcEqDistsThermal<nDim, nDist>(t, velocity, eqDist.data()); }
  IF_CONSTEXPR(thermalMode == 1) { lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(t, rho, velocity, eqDist.data()); }
  IF_CONSTEXPR(thermalMode == 2) { lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(t, rho, velocity, eqDist.data()); }
#endif
  return eqDist;
}

/**
 * \brief Calls function to return the thermal equilibrium distribution
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Moritz Waldmann
 *
 * \param[in] t                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDistsThermal(const MFloat t, const MFloat rho,
                                                                               const MFloat squaredVelocity,
                                                                               const MFloat* const velocity) {
  if(m_innerEnergy) {
    return getEqDistsThermal_<1>(t, rho, squaredVelocity, velocity);
  } else if(m_totalEnergy) {
    return getEqDistsThermal_<2>(t, rho, squaredVelocity, velocity);
  } else {
    return getEqDistsThermal_<0>(t, rho, squaredVelocity, velocity);
  }
}

/**
 * \brief Return thermal distributions to equilibrium
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Moritz Waldmann
 *
 * \tparam thermalMode Determines which equilibrium formulation is used
 *
 * \param[in] t                 Macroscopic temperature
 * \param[in] rho               Macroscopic density
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt thermalMode>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDistsThermal_(const MFloat t, const MFloat rho,
                                                                                const MFloat squaredVelocity,
                                                                                const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  IF_CONSTEXPR(thermalMode == 0) {
    lbfunc::calcEqDistsThermal<nDim, nDist>(t, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                            m_tp.data(), m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 1) {
    lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(t, rho, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(),
                                                m_mFld2.data(), m_tp.data(), m_distFld.data());
  }
  IF_CONSTEXPR(thermalMode == 2) {
    lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(t, rho, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(),
                                                m_mFld2.data(), m_tp.data(), m_distFld.data());
  }
#else
  IF_CONSTEXPR(thermalMode == 0) {
    lbfunc::calcEqDistsThermal<nDim, nDist>(t, squaredVelocity, velocity, eqDist.data());
  }
  IF_CONSTEXPR(thermalMode == 1) {
    lbfunc::calcEqDistsInnerEnergy<nDim, nDist>(t, rho, squaredVelocity, velocity, eqDist.data());
  }
  IF_CONSTEXPR(thermalMode == 2) {
    lbfunc::calcEqDistsTotalEnergy<nDim, nDist>(t, rho, squaredVelocity, velocity, eqDist.data());
  }
#endif
  return eqDist;
}

/**
 * \brief Return transport distributions to equilibrium
 *
 * Use this version if the squaredVelocity is not already precomputed
 *
 * \author Moritz Waldmann
 *
 * \param[in] c                 Macroscopic transport quantity
 * \param[in] rho               Macroscopic density
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDistsTransport(const MFloat c,
                                                                                 const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  lbfunc::calcEqDistsTransport<nDim, nDist>(c, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                                            m_distFld.data());
#else
  lbfunc::calcEqDistsTransport<nDim, nDist>(c, velocity, eqDist.data());
#endif
  return eqDist;
}

/**
 * \brief Return transport distributions to equilibrium
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Moritz Waldmann
 *
 * \param[in] c                 Macroscopic transport quantity
 * \param[in] rho               Macroscopic density
 * \param[in] squaredVelocity   Macroscopic velocity squared
 * \param[in] velocity          Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
std::array<MFloat, nDist> LbSolverDxQy<nDim, nDist, SysEqn>::getEqDistsTransport(const MFloat c,
                                                                                 const MFloat squaredVelocity,
                                                                                 const MFloat* const velocity) {
  std::array<MFloat, nDist> eqDist;
#ifdef WAR_NVHPC_PSTL
  lbfunc::calcEqDistsTransport<nDim, nDist>(c, squaredVelocity, velocity, eqDist.data(), m_mFld1.data(), m_mFld2.data(),
                                            m_tp.data(), m_distFld.data());
#else
  lbfunc::calcEqDistsTransport<nDim, nDist>(c, squaredVelocity, velocity, eqDist.data());
#endif
  return eqDist;
}

/**
 * \brief Calculate macroscopic variables for a given cell
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  cellId Cell id for which the macroscopic variables are calculated
 * \param[out] rho    Macroscopic density
 * \param[out] u      Macroscopic velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool compressible>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::calculateMacroscopicVariables(const MInt cellId, MFloat& rho,
                                                                             MFloat* const u) {
#ifdef WAR_NVHPC_PSTL
  std::array<MFloat, nDist> oldDist;
  for(MInt dist = 0; dist < nDist; dist++) {
    oldDist[dist] = a_oldDistribution(cellId, dist);
  }
  const MInt fldlen = Ld::dxQyFld();
  lbfunc::calcMacroVars<nDim, nDist, compressible>(oldDist.data(), rho, u, m_pFld.data(), m_nFld.data(), fldlen);
#else
  const MFloat* const dist = &a_oldDistribution(cellId, 0);
  lbfunc::calcMacroVars<nDim, nDist, compressible>(dist, rho, u);
#endif
}

/** \brief Calculate and set momentum flux for a given cell
 *  \author Miro Gondrum
 *  \param[in]  pCellId       Cell id for which the momentum flux is calculated
 *  \param[out] momentumFlux  momentumFlux tensor
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool compressible>
inline void LbSolverDxQy<nDim, nDist, SysEqn>::calculateMomentumFlux(const MInt pCellId, MFloat* const momentumFlux) {
  TRACE();

  const MFloat rho = a_variable(pCellId, PV->RHO);
#ifndef WAR_NVHPC_PSTL
  const MFloat* const u = &a_variable(pCellId, PV->U);
  const MFloat* const dist = &a_oldDistribution(pCellId, 0);
  sysEqn().calcMomentumFlux(rho, u, dist, momentumFlux);
#else
  std::array<MFloat, nDim> u;
  std::array<MFloat, nDist> dist;
  for(MInt d = 0; d < nDim; d++) {
    u[d] = a_variable(pCellId, d);
  }
  for(MInt d = 0; d < nDist; d++) {
    dist[d] = a_oldDistribution(pCellId, d);
  }
  sysEqn().calcMomentumFlux(rho, u.data(), dist.data(), momentumFlux);
#endif
}

/**
 * \brief Calculate and set momentum flux for a given cell
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] pCellId Cell id for which the momentum flux is calculated
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool compressible>
void LbSolverDxQy<nDim, nDist, SysEqn>::calculateMomentumFlux(const MInt pCellId) {
  TRACE();

  calculateMomentumFlux<compressible>(pCellId, m_momentumFlux[pCellId]);
}

#endif
