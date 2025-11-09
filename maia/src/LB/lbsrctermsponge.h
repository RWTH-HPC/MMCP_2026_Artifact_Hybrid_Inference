// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <unordered_map>
#include "COMM/mpioverride.h"
#include "IO/context.h"
#include "UTIL/debug.h"
#include "lbsolverdxqy.h"
#include "lbsrcterm.h"
#include "lbsrctermcontroller.h"

namespace maia::lb {

/** \brief  Base class for lb sponge layers
 *  \author Miro Gondrum
 *  \date   22.12.2023
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTerm_sponge : public LbSrcTerm<nDim, nDist, SysEqn> {
 private:
  void updateMapCellId2PmlCellId();

 protected:
  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;

  MBool m_isActive = false;

  struct SpongeCell {
    MInt cellId = -1;
    MFloat spongeFactor = -1.0;
  };

  std::vector<SpongeCell> m_spongeCells;
  std::unordered_map<MInt, MInt> m_mapCellId2SpongeCellId;

  MFloat m_spongeLayerThickness = 0.0;
  std::array<MFloat, 2 * nDim> m_spongeThicknessFactor{};
  MFloat m_spongeSigma = 0.0;
  MFloat m_spongeBeta = 2.0;
  MFloat m_trgRho = 1.0;
  std::array<MFloat, nDim> m_trgU{};

  /** \brief  Get index of cellId in m_spongeCells, returns -1 if unavailable
   *  \author Miro Gondrum
   *  \date   22.12.2023
   */
  MInt a_spongeCellId(const MInt cellId) {
    MInt spongeCellId = -1;
    auto search = m_mapCellId2SpongeCellId.find(cellId);
    if(search != m_mapCellId2SpongeCellId.end()) {
      spongeCellId = search->second;
    }
    return spongeCellId;
  };

  void readProperties() override;

 public:
  LbSrcTerm_sponge(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : m_solver(p_solver) { readProperties(); };

  void init() override;
};

/** \brief  Update m_mapCellIds2PmlCellId
 *  \author Miro Gondrum
 *  \date   22.12.2023
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_sponge<nDim, nDist, SysEqn>::updateMapCellId2PmlCellId() {
  const MInt noSpongeCells = m_spongeCells.size();
  m_mapCellId2SpongeCellId.clear();
  m_mapCellId2SpongeCellId.reserve(noSpongeCells);
  for(MInt spongeCellId = 0; spongeCellId < noSpongeCells; spongeCellId++) {
    m_mapCellId2SpongeCellId.emplace(m_spongeCells[spongeCellId].cellId, spongeCellId);
  }
}

/** \brief  Reading basics properties common for all lb sponge terms
 *  \author Miro Gondrum
 *  \date   22.12.2023
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_sponge<nDim, nDist, SysEqn>::readProperties() {
  TRACE();
  const MInt solverId = m_solver->m_solverId;

  std::stringstream ss;
  ss << "--------------------------------------------------------------------------------" << std::endl;
  ss << " Sponge properties" << std::endl;
  ss << "--------------------------------------------------------------------------------" << std::endl;
  ss << " [Sponge region]" << std::endl;

  // sponge information
  m_spongeSigma = Context::getSolverProperty<MFloat>("spongeSigma", solverId, AT_, &m_spongeSigma);
  m_spongeLayerThickness =
      Context::getSolverProperty<MFloat>("spongeLayerThickness", solverId, AT_, &m_spongeLayerThickness);
  ss << "  spongeSigma:             " << m_spongeSigma << std::endl;
  ss << "  spongeLayerThickness:    " << m_spongeLayerThickness << std::endl;
  ss << "  spongeThicknessFactor: [ ";
  for(MInt i = 0; i < 2 * nDim; i++) {
    m_spongeThicknessFactor[i] = Context::getSolverProperty<MFloat>("spongeThicknessFactor", solverId, AT_, i);
    ss << m_spongeThicknessFactor[i] << " ";
  }
  ss << "]" << std::endl;
  ss << " [Sponge shape]" << std::endl;
  m_spongeBeta = Context::getSolverProperty<MFloat>("spongeBeta", solverId, AT_, &m_spongeBeta);
  ss << "  spongeBeta:              " << m_spongeBeta << std::endl;
  ss << " [Sponge target state]" << std::endl;
  m_trgRho = (m_solver->m_densityFluctuations) ? 0.0 : 1.0;
  m_trgRho = Context::getSolverProperty<MFloat>("spongeTrgRho", solverId, AT_, &m_trgRho);
  ss << "  spongeTrgRho:            " << m_trgRho << std::endl;
  m_trgU.fill(0.0);
  for(MInt i = 0; i < nDim; i++) {
    m_trgU[i] = Context::getSolverProperty<MFloat>("spongeTrgU", solverId, AT_, &m_trgU[i], i) * m_solver->m_Ma * LBCS;
  }
  ss << "  spongeTrgU:              ";
  for(MInt i = 0; i < nDim; i++) {
    ss << m_trgU[i] << " ";
  }
  ss << std::endl;
  ss << "--------------------------------------------------------------------------------" << std::endl;

  if(m_solver->domainId() == 0) std::cout << ss.str();
  m_log << ss.str();
}

/** \brief  Initialize properties common by all lb sponge terms
 *  \author Miro Gondrum
 *  \date   22.12.2023
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_sponge<nDim, nDist, SysEqn>::init() {
  TRACE();
  m_isActive = (m_spongeLayerThickness > 0.0) && m_solver->isActive();
  if(!m_isActive) return;

  const MFloat epsilon = m_solver->c_cellLengthAtLevel(m_solver->maxLevel()) / 1000.0;
  std::array<MFloat, 2 * nDim> domainBoundaries{};
  {
    const MFloat halfCellWidth = F1B2 * m_solver->grid().cellLengthAtCell(0);
    for(MInt i = 0; i < nDim; i++) {
      domainBoundaries[2 * i + 0] = m_solver->a_coordinate(0, i) - halfCellWidth;
      domainBoundaries[2 * i + 1] = m_solver->a_coordinate(0, i) + halfCellWidth;
    }
  }

  for(MInt cellId = 0; cellId < m_solver->m_cells.size(); cellId++) {
    if(!m_solver->c_isLeafCell(cellId)) continue;
    const MFloat halfCellWidth = F1B2 * m_solver->grid().cellLengthAtCell(cellId);
    for(MInt i = 0; i < nDim; i++) {
      domainBoundaries[2 * i + 0] =
          std::min(domainBoundaries[2 * i + 0], m_solver->a_coordinate(cellId, i) - halfCellWidth);
      domainBoundaries[2 * i + 1] =
          std::max(domainBoundaries[2 * i + 1], m_solver->a_coordinate(cellId, i) + halfCellWidth);
    }
  }
  // exchange domain boundaries
  std::array<MFloat, nDim> tmpMin{};
  std::array<MFloat, nDim> tmpMax{};
  for(MInt i = 0; i < nDim; i++) {
    tmpMin[i] = domainBoundaries[2 * i + 0];
    tmpMax[i] = domainBoundaries[2 * i + 1];
  }
  MPI_Allreduce(MPI_IN_PLACE, tmpMin.data(), nDim, maia::type_traits<MFloat>::mpiType(), MPI_MIN, m_solver->mpiComm(),
                AT_, "MPI_IN_PLACE", "tmpMin");
  MPI_Allreduce(MPI_IN_PLACE, tmpMax.data(), nDim, maia::type_traits<MFloat>::mpiType(), MPI_MAX, m_solver->mpiComm(),
                AT_, "MPI_IN_PLACE", "tmpMax");
  for(MInt i = 0; i < nDim; i++) {
    domainBoundaries[2 * i + 0] = tmpMin[i];
    domainBoundaries[2 * i + 1] = tmpMax[i];
  }

  // fill sponge coordinate
  std::array<MFloat, 2 * nDim> spongeCoord{};
  for(MInt i = 0; i < nDim; i++) {
    spongeCoord[2 * i + 0] = domainBoundaries[2 * i + 0] + m_spongeLayerThickness * m_spongeThicknessFactor[2 * i + 0];
    spongeCoord[2 * i + 1] = domainBoundaries[2 * i + 1] - m_spongeLayerThickness * m_spongeThicknessFactor[2 * i + 1];
  }
  std::array<MFloat, 2 * nDim> F1BspongeCoord{};
  for(MInt i = 0; i < nDim; i++) {
    F1BspongeCoord[2 * i + 0] = 1.0 / std::max(epsilon, m_spongeLayerThickness * m_spongeThicknessFactor[2 * i + 0]);
    F1BspongeCoord[2 * i + 1] = 1.0 / std::max(epsilon, m_spongeLayerThickness * m_spongeThicknessFactor[2 * i + 1]);
  }

  // search sponge ids
  for(MInt cellId = 0; cellId < m_solver->m_cells.size(); cellId++) {
    if(m_solver->a_isHalo(cellId)) continue;

    std::array<MFloat, nDim> dist{};
    for(MInt i = 0; i < nDim; i++) {
      const MUint distIdNegSide = 2 * i + 0;
      const MUint distIdPosSide = 2 * i + 1;
      const MFloat distNegSide =
          (spongeCoord[distIdNegSide] - m_solver->a_coordinate(cellId, i)) * F1BspongeCoord[distIdNegSide];
      const MFloat distPosSide =
          (m_solver->a_coordinate(cellId, i) - spongeCoord[distIdPosSide]) * F1BspongeCoord[distIdPosSide];
      dist[i] = std::max(distNegSide, distPosSide);
    }

    // choose the farthest sponge region
    const MFloat relDist = *std::max_element(dist.begin(), dist.end());

    if(relDist <= 0.0) continue;

    const MFloat nodeDist = std::pow(relDist, m_spongeBeta);
    const MFloat spongeFactor = m_spongeSigma * nodeDist;

    constexpr MFloat epsSpongeFactor = 1e-15;
    if(spongeFactor > epsSpongeFactor) {
      auto& cell = m_spongeCells.emplace_back();
      cell.cellId = cellId;
      cell.spongeFactor = spongeFactor;
    }
  }
  // Update map cellId->spongeCellId
  updateMapCellId2PmlCellId();
}

} // namespace maia::lb
