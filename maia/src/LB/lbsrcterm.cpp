// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsrcterm.h"

#include "IO/context.h"
#include "UTIL/debug.h"
#include "lbsolverdxqy.h"
#include "lbsrctermcontroller.h"

#include <vector>

namespace maia::lb {

//---LbSrcTerm_monopole-----------------------------------------------------
//-declaration
/** \brief  Class to handle acoustic monopole source terms in lb solver
 *  \author Miro Gondrum, Benyamin Krisna
 *  \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTerm_monopole : public LbSrcTerm<nDim, nDist, SysEqn> {
 private:
  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;

  struct MonopoleInfo {
    MInt id = -1;
    MBool isActive = false;
    MFloat amplitude; // of pulsating mass
    MFloat rhoFluct;
    MFloat radius = -1.0;
    MFloat strouhal;
    MFloat omega;
    MFloat phaseShift = 0.0;
    MBool windowing = false;
    std::array<MFloat, nDim> position{};
    std::vector<MInt> cellIds{};
    std::vector<MFloat> srcTerms{};
  };
  std::vector<MonopoleInfo> m_monopole;

 protected:
  void readProperties() override;

 public:
  //  static const LbRegSrcTerm<nDim, nDist, LbSrcTerm_monopole<nDim, nDist>> reg;

  LbSrcTerm_monopole(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : m_solver(p_solver) { readProperties(); };

  void init() override;
  void apply_preCollision() override{/*UNUSED FOR FIRST ORDER IMPLEMENTATION*/};
  void apply_postCollision() override;
  void apply_postPropagation() override{};
};

//-definiton
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_monopole<nDim, nDist, SysEqn>::readProperties() {
  TRACE();
  const MInt solverId = m_solver->m_solverId;
  MInt noMonopole = 0;
  for(MInt id = 0;; id++) {
    if(!Context::propertyExists("monopoleAmplitude_" + std::to_string(id))) {
      break;
    }
    noMonopole++;
    auto& monopole = m_monopole.emplace_back();
    monopole.id = id;
    for(MInt d = 0; d < nDim; d++) {
      monopole.position[d] = Context::getSolverProperty<MFloat>("monopolePos_" + std::to_string(id), solverId, AT_);
    }
    monopole.amplitude = Context::getSolverProperty<MFloat>("monopoleAmplitude_" + std::to_string(id), solverId, AT_);
    monopole.strouhal = Context::getSolverProperty<MFloat>("monopoleStrouhal_" + std::to_string(id), solverId, AT_);
    monopole.radius =
        Context::getSolverProperty<MFloat>("monopoleRadius_" + std::to_string(id), solverId, AT_, &monopole.radius);
    monopole.windowing = Context::getSolverProperty<MBool>("monopoleWindowing_" + std::to_string(id), solverId, AT_,
                                                           &monopole.windowing);
    const MFloat phaseShiftStrouhal = Context::getSolverProperty<MFloat>("monopolePhaseShift_" + std::to_string(id),
                                                                         solverId, AT_, &monopole.phaseShift);
    if(monopole.radius < m_solver->c_cellLengthAtLevel(m_solver->maxLevel())) {
      m_log << "monopoleRadius < max level cell length: Hence, it is clipped to max level cell length" << std::endl;
      monopole.radius = m_solver->c_cellLengthAtLevel(m_solver->maxLevel());
    }
    const MFloat u_infty = m_solver->m_Ma * LBCS;
    const MFloat dxLb = m_solver->c_cellLengthAtLevel(m_solver->maxLevel());
    const MFloat massAmp = monopole.amplitude * (POW3(1.0 / dxLb)); // from STL to LB units m=rho*L^3
    const MFloat densityAmp = massAmp / (F4B3 * PI * POW3(monopole.radius / dxLb));
    monopole.omega = 2.0 * PI * (monopole.strouhal * u_infty / m_solver->m_referenceLength);
    monopole.rhoFluct = densityAmp;
    monopole.phaseShift = 2.0 * PI * (phaseShiftStrouhal * u_infty / m_solver->m_referenceLength);
    if(m_solver->domainId() == 0) {
      std::stringstream ss;
      ss << "INFO: LbSrcTerm_monopole " << monopole.id << " analytical solution:" << std::endl;
      ss << "rho(r,t)= 1.0 -" << (massAmp * POW2(monopole.omega)) / (4 * PI * POW2(LBCS)) << "/ (mag(coords)/" << dxLb
         << ") * cos(" << monopole.omega << "* (" << g_timeSteps - 1 << " - mag(coords)/" << dxLb << "/" << LBCS << "))"
         << std::endl;
      std::cout << ss.str();
    }
  }
  if(noMonopole == 0) {
    TERMM(1, "Property monopoleAmplitude_0 is missing!");
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_monopole<nDim, nDist, SysEqn>::init() {
  TRACE();
  for(auto& monopole : m_monopole) {
    // First, clear everything -> usable to reinit
    monopole.cellIds.clear();
    monopole.srcTerms.clear();
    // Now initialize
    if(monopole.radius < 0.0) {
      const MInt monopoleCellId = m_solver->getIdAtPoint(monopole.position.data());
      if(monopoleCellId != -1) {
        monopole.isActive = true;
        monopole.cellIds.push_back(monopoleCellId);
      }
    } else {
      const MFloat radiusSq = POW2(monopole.radius);
      for(MInt i = 0; i < m_solver->m_currentMaxNoCells; i++) {
        const MInt cellId = m_solver->m_activeCellList[i];
        MFloat radiusSq_ = 0.0;
        for(MInt d = 0; d < nDim; d++) {
          radiusSq_ += POW2(monopole.position[d] - m_solver->a_coordinate(cellId, d));
        }
        if(radiusSq_ < radiusSq) {
          monopole.isActive = true;
          monopole.cellIds.push_back(cellId);
        }
      }
    }
    const MInt noSrcTermCells = monopole.cellIds.size();
    monopole.srcTerms.resize(noSrcTermCells, 0.0);
    //- info output
    MInt noSrcTermCellsGlobal = noSrcTermCells;
    MPI_Allreduce(MPI_IN_PLACE, &noSrcTermCellsGlobal, 1, MPI_INT, MPI_SUM, m_solver->mpiComm(), AT_, "MPI_IN_PLACE",
                  "noSrcTermCellsGlobal");
    const MFloat ppw = 2 * PI * LBCS / monopole.omega; // ppw = lambda in LB units
    const MFloat periode = m_solver->m_referenceLength / (monopole.strouhal * m_solver->m_Ma * LBCS);
    std::stringstream ss;
    ss << "INFO: LbSrcTerm_monopole " << monopole.id << ":" << std::endl;
    ss << "    amplitude             :" << monopole.amplitude << std::endl;
    ss << "    rhoFluct              :" << monopole.rhoFluct << std::endl;
    ss << "    radius                :" << monopole.radius << std::endl;
    ss << "    strouhal              :" << monopole.strouhal << std::endl;
    ss << "    omega                 :" << monopole.omega << std::endl;
    ss << "    phaseShift            :" << monopole.phaseShift << std::endl;
    ss << "    windowing             :" << monopole.windowing << std::endl;
    ss << "    noSrcTermCells        :" << noSrcTermCellsGlobal << std::endl;
    ss << "    points per wave       :" << ppw << std::endl;
    ss << "    timesteps per periode :" << periode << std::endl;
    if(m_solver->domainId() == 0) {
      std::cout << ss.str();
    }
    m_log << ss.str();
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_monopole<nDim, nDist, SysEqn>::apply_postCollision() {
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  TRACE();
  for(auto& monopole : m_monopole) {
    if(!monopole.isActive) continue;
    const MInt noSrcTermCells = monopole.cellIds.size();
    for(MInt i = 0; i < noSrcTermCells; i++) {
      //-update source term
      monopole.srcTerms[i] =
          -monopole.rhoFluct * monopole.omega * sin((monopole.omega + monopole.phaseShift) * globalTimeStep);
      if(monopole.windowing) {
        // Windowing shall be performed for the first wave length
        if(globalTimeStep < 2 * PI / (monopole.omega)) {
          const MFloat S = (globalTimeStep == 0 ? 0.5 : 1);
          const MFloat hannHalfWindow = S * 0.5 * (1 - cos(monopole.omega * globalTimeStep / 2.0));
          monopole.srcTerms[i] *= hannHalfWindow;
        }
      }

      //-apply source term in solver
      const MInt cellId = monopole.cellIds[i];
      for(MInt j = 0; j < nDist; j++) {
        m_solver->a_distribution(cellId, j) += Ld::tp(Ld::distType(j)) * monopole.srcTerms[i];
      }
    }
  }
}

//-registration for LbSrcTermFactory class
// V1 : no  static member reg in LbSrcTerm_monopole needed
static const LbRegSrcTerm<LbSrcTerm_monopole<2, 9, LbSysEqnIncompressible<2, 9>>> reg_Monopoled2q9("LB_SRC_MONOPOLE");
static const LbRegSrcTerm<LbSrcTerm_monopole<3, 19, LbSysEqnIncompressible<3, 19>>>
    reg_Monopoled3q19("LB_SRC_MONOPOLE");
static const LbRegSrcTerm<LbSrcTerm_monopole<3, 27, LbSysEqnIncompressible<3, 27>>>
    reg_Monopoled3q27("LB_SRC_MONOPOLE");
static const LbRegSrcTerm<LbSrcTerm_monopole<2, 9, LbSysEqnCompressible<2, 9>>> reg_Monopoled2q9C("LB_SRC_MONOPOLE");
static const LbRegSrcTerm<LbSrcTerm_monopole<3, 19, LbSysEqnCompressible<3, 19>>> reg_Monopoled3q19C("LB_SRC_MONOPOLE");
static const LbRegSrcTerm<LbSrcTerm_monopole<3, 27, LbSysEqnCompressible<3, 27>>> reg_Monopoled3q27C("LB_SRC_MONOPOLE");
// V2 : using  static member reg in LbSrcTerm_monopole
// template <MInt nDim, MInt nDist, class SysEqn>
// LbRegSrcTerm<nDim, nDist, LbSrcTerm_monopole<nDim, nDist>>
//    LbSrcTerm_monopole<nDim, nDist>::reg("LB_SRC_MONOPOLE");
// template class LbSrcTerm_monopole<3, 19>;
// template class LbSrcTerm_monopole<3, 27>;

//---LbSrcTerm_XXX----------------------------------------------------------
// add furhter source terms below

} // namespace maia::lb
