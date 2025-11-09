// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsrctermsponge.h"

#include "UTIL/parallelfor.h"
#include "lbsrctermcontroller.h"
#include "lbsyseqn.h"

namespace maia::lb {

//---LbSrcTerm_spongeRhoConst---------------------------------------------------
/** \brief LB sponge source term : sponge towards const target f_eq(rho_trg, u)
 *  \author Miro Gondrum
 *  \date   04.01.2024
 *  The PPDF f is damped towards an equilibrium state f_a following:
 *    f -= sigma *(f_eq - f_a) , with f_a = f_eq(rho_trg, u)
 *  Here, rho_trg is a user defined input and for this sponge are constant in
 *  space and time and u is the local velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTerm_spongeRhoConst : public LbSrcTerm_sponge<nDim, nDist, SysEqn> {
 public:
  using Base = LbSrcTerm_sponge<nDim, nDist, SysEqn>;
  LbSrcTerm_spongeRhoConst(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : Base(p_solver){};
  void apply_preCollision() override{};
  void apply_postCollision() override;
  void apply_postPropagation() override{};

 private:
};

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_spongeRhoConst<nDim, nDist, SysEqn>::apply_postCollision() {
  TRACE();
  LbSolverDxQy<nDim, nDist, SysEqn>* const s = this->m_solver; // alias for readability
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  const MInt timeStep = s->getCurrentTimeStep();
  const MUint noSpongeCells = this->m_spongeCells.size();
  maia::parallelFor<true>(0, noSpongeCells, [=](MInt spongeCellId) {
    auto& cell = this->m_spongeCells[spongeCellId];
    const MInt cellId = cell.cellId;
    if((timeStep - 1) % IPOW2(s->maxLevel() - s->a_level(cellId)) != 0) return;
    // compute the forcing term for the density
    const MFloat deltaRho = s->a_variable(cellId, s->PV->RHO) - this->m_trgRho;
    s->a_variable(cellId, s->PV->RHO) -= cell.spongeFactor * deltaRho;

    // prescribe new density
    s->a_distribution(cellId, Ld::lastId()) = s->a_variable(cellId, s->PV->RHO);
    for(MInt dist = 0; dist < nDist - 1; dist++) {
      s->a_distribution(cellId, nDist - 1) -= s->a_distribution(cellId, dist);
    }
  });
}

//---LbSrcTerm_spongeEqConst----------------------------------------------------
/** \brief LB sponge source term : sponge towards const target f_eq(rho_trg, u_trg)
 *  \author Miro Gondrum
 *  \date   04.01.2024
 *  The PPDF f is damped towards a predefined equilibrium state f_a following:
 *    f -= sigma *(f_eq - f_a) , with f_a = f_eq(rho_trg, u_trg)
 *  Here, rho_trg and u_trg are user defined input and for this sponge are
 *  constant in space and time.
 *
 *  Reference:
 *    - Kam et al. 2007: http://dx.doi.org/10.2514/1.27632 (Type 3: ABC)
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTerm_spongeEqConst : public LbSrcTerm_sponge<nDim, nDist, SysEqn> {
 public:
  using Base = LbSrcTerm_sponge<nDim, nDist, SysEqn>;
  LbSrcTerm_spongeEqConst(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : Base(p_solver){};
  void init() override;
  void apply_preCollision() override{};
  void apply_postCollision() override;
  void apply_postPropagation() override{};

 private:
  std::array<MFloat, nDist> m_trgEq{};
};

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_spongeEqConst<nDim, nDist, SysEqn>::init() {
  LbSolverDxQy<nDim, nDist, SysEqn>& s = *(this->m_solver); // alias for readability
  Base::init();
  m_trgEq = s.getEqDists(this->m_trgRho, this->m_trgU.data());
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_spongeEqConst<nDim, nDist, SysEqn>::apply_postCollision() {
  TRACE();
  LbSolverDxQy<nDim, nDist, SysEqn>* const s = this->m_solver; // alias for readability

  const MInt timeStep = s->getCurrentTimeStep();
  const MUint noSpongeCells = this->m_spongeCells.size();
  maia::parallelFor<true>(0, noSpongeCells, [=](MInt spongeCellId) {
    auto& cell = this->m_spongeCells[spongeCellId];
    const MInt cellId = cell.cellId;
    if((timeStep - 1) % IPOW2(s->maxLevel() - s->a_level(cellId)) != 0) return;
    std::array<MFloat, nDim> u{};
    for(MInt d = 0; d < nDim; d++) {
      u[d] = s->a_variable(cellId, d);
    }
    std::array<MFloat, nDist> eqDist = s->getEqDists(s->a_variable(cellId, s->PV->RHO), u.data());
    for(MInt i = 0; i < nDist; i++) {
      s->a_distribution(cellId, i) -= cell.spongeFactor * (eqDist[i] - m_trgEq[i]);
    }
  });
}

//---LbSrcTerm_spongeViscosity--------------------------------------------------
/** \brief LB sponge source term : sponge viscosity towards outside
 *  \author Miro Gondrum
 *  \date   04.01.2024
 *  The PPDF f is damped towards a its equilibrium state f_eq following:
 *    f -= omega_sponge *(f_eq - f)
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTerm_spongeVisocity : public LbSrcTerm_sponge<nDim, nDist, SysEqn> {
 public:
  using Base = LbSrcTerm_sponge<nDim, nDist, SysEqn>;
  LbSrcTerm_spongeVisocity(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : Base(p_solver){};
  void apply_preCollision() override;
  void apply_postCollision() override{};
  void apply_postPropagation() override{};

 private:
};

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_spongeVisocity<nDim, nDist, SysEqn>::apply_preCollision() {
  TRACE();
  LbSolverDxQy<nDim, nDist, SysEqn>* const s = this->m_solver; // alias for readability

  const MInt timeStep = s->getCurrentTimeStep();
  const MUint noSpongeCells = this->m_spongeCells.size();
  maia::parallelFor<true>(0, noSpongeCells, [=](MInt spongeCellId) {
    auto& cell = this->m_spongeCells[spongeCellId];
    const MInt cellId = cell.cellId;
    if((timeStep - 1) % IPOW2(s->maxLevel() - s->a_level(cellId)) != 0) return;
    s->a_nu(cellId) *= (1.0 + cell.spongeFactor);
  });
}

//---registration for LbSrcTermFactory class------------------------------------
static const LbRegSrcTerm<LbSrcTerm_spongeRhoConst<2, 9, LbSysEqnCompressible<2, 9>>>
    reg_SpongeRhoConstd2q9C("LB_SPONGE_RHOCONST");
static const LbRegSrcTerm<LbSrcTerm_spongeRhoConst<2, 9, LbSysEqnIncompressible<2, 9>>>
    reg_SpongeRhoConstd2q9("LB_SPONGE_RHOCONST");
static const LbRegSrcTerm<LbSrcTerm_spongeRhoConst<3, 27, LbSysEqnCompressible<3, 27>>>
    reg_SpongeRhoConstd3q27C("LB_SPONGE_RHOCONST");
static const LbRegSrcTerm<LbSrcTerm_spongeRhoConst<3, 27, LbSysEqnIncompressible<3, 27>>>
    reg_SpongeRhoConstd3q27("LB_SPONGE_RHOCONST");

static const LbRegSrcTerm<LbSrcTerm_spongeEqConst<2, 9, LbSysEqnCompressible<2, 9>>>
    reg_SpongeEqConstd2q9C("LB_SPONGE_EQCONST");
static const LbRegSrcTerm<LbSrcTerm_spongeEqConst<2, 9, LbSysEqnIncompressible<2, 9>>>
    reg_SpongeEqConstd2q9("LB_SPONGE_EQCONST");
static const LbRegSrcTerm<LbSrcTerm_spongeEqConst<3, 27, LbSysEqnCompressible<3, 27>>>
    reg_SpongeEqConstd3q27C("LB_SPONGE_EQCONST");
static const LbRegSrcTerm<LbSrcTerm_spongeEqConst<3, 27, LbSysEqnIncompressible<3, 27>>>
    reg_SpongeEqConstd3q27("LB_SPONGE_EQCONST");

static const LbRegSrcTerm<LbSrcTerm_spongeVisocity<2, 9, LbSysEqnCompressible<2, 9>>>
    reg_SpongeViscoisty2q9C("LB_SPONGE_VISCOSITY");
static const LbRegSrcTerm<LbSrcTerm_spongeVisocity<2, 9, LbSysEqnIncompressible<2, 9>>>
    reg_SpongeViscoisty2q9("LB_SPONGE_VISCOSITY");

} // namespace maia::lb
