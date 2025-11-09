// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lptbase.h"
#include "INCLUDE/maiaconstants.h"
#include "UTIL/functions.h"
#include "lpt.h"


template <MInt nDim>
MInt LPTBase<nDim>::s_interpolationOrder = 0;
template <MInt nDim>
MInt LPTBase<nDim>::s_interpolationMethod = 0;
template <MInt nDim>
MFloat LPTBase<nDim>::s_distFactorImp = 0.0;
template <MInt nDim>
LPT<nDim>* LPTBase<nDim>::s_backPtr = nullptr;

template <MInt nDim>
template <MInt a, MInt b>
void LPTBase<nDim>::interpolateAndCalcWeights(const MInt cellId, const MFloat* const x, MFloat* const result,
                                              std::vector<MFloat>& weight, MFloat* const gradientResult) {
  std::vector<MInt>* nghbrL = nullptr;
  std::vector<MInt> tmp;
  if(cellId == m_cellId) {
    getNghbrList(m_neighborList, cellId);
    nghbrL = &m_neighborList;
  } else {
    getNghbrList(tmp, cellId);
    nghbrL = &tmp;
  }


  // calculate weights
  if(s_interpolationOrder > 0 || s_backPtr->m_couplingRedist) {
    // calculate weights of the redistribution
    const auto noNghbr = static_cast<MUint>(nghbrL->size());

    // higher values means a sharper fall off
    // static constexpr MFloat distFactorImp = 0.5;
    // static constexpr MInt s_interpolationMethod = 2;//default was 1;

    if(noNghbr > 0 && (s_interpolationMethod > 0 || noNghbr != weight.size())) {
      MFloat sumImp = 0.0;
      weight.clear();
      for(MUint n = 0; n < noNghbr; n++) {
        MFloat dx = 0;
        if(s_interpolationMethod > 0) {
          for(MInt i = 0; i < nDim; i++) {
            const MFloat coord = s_backPtr->a_coordinate(nghbrL->at(n), i);
            dx += POW2(coord - x[i]) + MFloatEps;
          }
        }

        switch(s_interpolationMethod) {
          case 0:
            // average overall cells
            weight.push_back(1.0);
            break;
          case 1:
            // linear
            weight.push_back(1.0 / dx);
            break;
          case 2: {
            // normal distribution
            const MFloat cellLength = s_backPtr->c_cellLengthAtLevel(s_backPtr->c_level(m_cellId));
            // higher values means a sharper fall off
            MFloat distFactorImp = s_distFactorImp;
            if(s_backPtr->m_engineSetup && s_backPtr->a_coordinate(m_cellId, 1) > 0.05) {
              distFactorImp = 0.25;
            }
            weight.push_back(exp(-distFactorImp * dx / POW2(cellLength)));
            break;
          }
          default:
            TERMM(-1, "Invalid interpolation method.");
        }
        sumImp += weight.back();
        ASSERT(!std::isnan(weight.back()), "ERROR: Weight is NaN with dx = " + std::to_string(dx) + "!");
      }

      // normalize weights
      for(MUint n = 0; n < noNghbr; n++) {
        weight.at(n) /= sumImp;
      }
    }
  }

  if(b == 0) return;

  MInt c = b;
  if(b > s_backPtr->PV.noVars()) {
    c = b - 1;
  }

  switch(s_interpolationOrder) {
    case 0:
      // use current cell value
      for(MInt i = a; i < c; i++) {
        result[i] = s_backPtr->a_fluidVariable(cellId, i);
      }
      if(b > s_backPtr->PV.noVars()) {
        result[c] = s_backPtr->a_fluidSpecies(cellId);
      }
      break;
    case 1: {
      // linear interpolation
      MUint n = 0;
      for(auto nghbrId : *nghbrL) {
        for(MInt i = a; i < c; i++) {
          result[i] += weight[n] * s_backPtr->a_fluidVariable(nghbrId, i);
        }
        if(b > s_backPtr->PV.noVars()) {
          result[c] += weight[n] * s_backPtr->a_fluidSpecies(nghbrId);
        }
        n++;
      }
      break;
    }
    case 2:
    case 3:
      // is used by ellipsoids
      // least-square interpolation
      {
        if(b > s_backPtr->PV.noVars()) mTerm(1, AT_, "Not implemented yet!");
        if(gradientResult) {
          std::array<MFloat, b + nDim * nDim> tempResults;
          s_backPtr->template interpolateVariablesLS<a, b, true>(cellId, x, tempResults.data());
          for(MInt i = a; i < c; i++) {
            result[i - a] = tempResults[i];
          }
          for(MInt i = 0; i < nDim * nDim; i++) {
            gradientResult[i] = tempResults[c + i];
          }
        } else
          s_backPtr->template interpolateVariablesLS<a, b>(cellId, x, result);
      }
      break;
    default:
      TERMM(-1, "Invalid interpolation order.");
  }
}


template void LPTBase<2>::interpolateAndCalcWeights<0, 0>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<2>::interpolateAndCalcWeights<0, 2>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<2>::interpolateAndCalcWeights<0, 3>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<2>::interpolateAndCalcWeights<0, 5>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<2>::interpolateAndCalcWeights<0, 6>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<3>::interpolateAndCalcWeights<0, 0>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<3>::interpolateAndCalcWeights<0, 3>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<3>::interpolateAndCalcWeights<0, 4>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<3>::interpolateAndCalcWeights<0, 6>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);
template void LPTBase<3>::interpolateAndCalcWeights<0, 7>(const MInt cellId, const MFloat* const x,
                                                          MFloat* const result, std::vector<MFloat>& weight,
                                                          MFloat* const gradientResult);

/// Get the neighborlist from the Lagrange particle container object
/// \param neighborList List of the neighboring cells
/// \param cellId Get neighbors of this cell
template <MInt nDim>
void LPTBase<nDim>::getNghbrList(std::vector<MInt>& neighborList, const MInt cellId) {
  ASSERT(cellId >= 0, "Invalid cellId!");

  static const MInt noOfRedistLayers = (s_backPtr->m_noRedistLayer > 0) ? s_backPtr->m_noRedistLayer : 2;

  MUint additionalLayers = 0;

  if(s_backPtr->a_volumeFraction(m_cellId) > 0.75) {
    // use further redistribution of source-terms for packed cells
    // TODO labels:LPT otherwise include particel-particle interaction and shift particles into neighboring cells!
    additionalLayers = 1;
  } else if(!m_neighborList.empty() && m_neighborList[0] == cellId && !s_backPtr->a_isBndryCell(cellId)) {
    // still current no need to update
    return;
  }


#ifdef _OPENMP
#pragma omp critical
#endif
  {
    if(s_backPtr->m_cellToNghbrHood.count(cellId) > 0
       && s_backPtr->m_cellToNghbrHood[cellId].size() >= 50 * additionalLayers) {
      // cell neighbor cells are already known, including additionalLayers
      neighborList = s_backPtr->m_cellToNghbrHood[cellId];
    } else {
      if(s_backPtr->m_cellToNghbrHood.count(cellId) > 0 && additionalLayers > 0) {
        // cell neighbors known, but increased to consider additional layers
        s_backPtr->m_cellToNghbrHood.erase(cellId);
      }
      neighborList.clear();

      // cells neighbor's are unknown
      std::vector<MInt> neighborListGrid;
      s_backPtr->grid().findNeighborHood(cellId, noOfRedistLayers + additionalLayers, neighborListGrid);
      for(MInt i = 0; i < (signed)neighborListGrid.size(); i++) {
        if(s_backPtr->a_isValidCell(neighborListGrid[i])) {
          neighborList.push_back(neighborListGrid[i]);
        }
      }
      s_backPtr->m_cellToNghbrHood.emplace(make_pair(cellId, neighborList));
    }
  }
}


/// \brief  Checks whether position is within cell with number cellId if not, cellId is updated.
///
/// \author Rudie Kunnen, Sven Berger, Tim Wegmann
template <MInt nDim>
void LPTBase<nDim>::checkCellChange(const MFloat* oldPosition, const MBool allowHaloNonLeaf) {
  TRACE();

  // no cell change when the position is identical
  if(oldPosition != nullptr
     && std::abs(m_position[0] - oldPosition[0]) + std::abs(m_position[1] - oldPosition[1])
                + std::abs(m_position[2] - oldPosition[2])
            < 1E-12) {
    return;
  }

  if(m_cellId != -1) m_cellId = s_backPtr->grid().findContainingLeafCell(m_position.data(), m_cellId, true);

  if(allowHaloNonLeaf && m_cellId < 0) {
    mTerm(1, AT_, "Invalid cell during injection! Send not possible, broadcast required!");
  }

  if(m_cellId > 0 && !s_backPtr->c_isLeafCell(m_cellId)) {
    if(!allowHaloNonLeaf) {
      m_cellId = -1;
    } else {
      reqSend() = true;
      isInvalid() = false;
      return;
    }
  }

  // a ghost particle may never become a real particle!
  auto partWasHalo = wasSend();
  // no valid cell found left geometry
  if(m_cellId < 0 || (!partWasHalo && !s_backPtr->a_isValidCell(m_cellId))) {
    toBeDeleted() = true;
    toBeRespawn() = true;
    isInvalid() = true;
    return;
  }

  // for multisolver, check whether the current cellId is a solver interface cell
  const MBool origCellIsHaloCell = s_backPtr->a_isHalo(m_cellId);
  const MBool origCellIsPeriodicHaloCell = s_backPtr->grid().isPeriodic(m_cellId);
  if(globalNoDomains() > 1) {
    // active particles may move to other solvers
    if(!isInvalid()) {
      // check whether particle changes solvers
      if(origCellIsHaloCell) {
        reqSend() = true;
        isWindow() = false;
      }

      // check whether particle changes solvers in periodic direction
      if(origCellIsPeriodicHaloCell) {
        isWindow() = false;
        isInvalid() = true;
        reqSend() = true;
      }

      isWindow() = s_backPtr->a_isWindow(m_cellId);
    }
  }
}


/// \brief Update particle properties
/// \author Tim Wegmann
/// \date   Aug 2021
template <MInt nDim>
void LPTBase<nDim>::updateProperties(const MBool init) {
  if(init) {
    initProperties();
  }


  if(s_backPtr->noDomains() > 1) {
    if(s_backPtr->a_isHalo(m_cellId)) {
      reqSend() = true;
      isInvalid() = true;
    } else if(s_backPtr->a_isWindow(m_cellId)) {
      isWindow() = true;
    }
  }
}

/// \brief particle-wall collision step
/// \author Sven Berger , Tim Wegmann, Laurent Andre
/// \date   June 2016, Jan 2021, Nov 2023
template <MInt nDim>
void LPTBase<nDim>::particleWallCollision() {
#ifdef LPT_DEBUG
  MLong debugPartId = -1;
#endif

  const MFloat dt = s_backPtr->m_timeStep - m_creationTime;

  const MInt bndId = s_backPtr->a_bndryCellId(m_cellId) > -1 ? s_backPtr->a_bndryCellId(m_cellId)
                                                             : s_backPtr->a_bndryCellId(m_oldCellId);

  if(bndId < 0) return;

  const MInt lastValidCellId =
      s_backPtr->a_isValidCell(m_cellId) ? m_cellId : s_backPtr->a_isValidCell(m_oldCellId) ? m_oldCellId : -1;

  const MFloat radius = effectiveWallCollisionRadius();

  std::array<MFloat, nDim> bodyVel{};
  MFloat nBodyVel = 0;
  for(MInt i = 0; i < nDim; i++) {
    nBodyVel += s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_velocity[i]
                * s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i];
  }
  if(nBodyVel > 0) {
    for(MInt i = 0; i < nDim; i++) {
      bodyVel[i] = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_velocity[i];
    }
  } else {
    for(MInt i = 0; i < nDim; i++) {
      bodyVel[i] = 0.0;
    }
  }

  MFloat nx0 = 0.0;
  MFloat na = 0.0;
  MFloat nb = 0.0;

  ASSERT(s_backPtr->m_bndryCells->a[bndId].m_noSrfcs == 1, "Not implemented for multiple surfaces yet!");

  std::array<MFloat, nDim> deltaX{};
  for(MInt i = 0; i < nDim; i++) {
    deltaX[i] = m_position[i] - m_oldPos[i];
    nx0 += s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i]
           * s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_coordinates[i];
    na += s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i] * m_oldPos[i];
    nb += s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i] * deltaX[i];
  }

  // We have a geometrical overlap if (x(t) - x_surf)*normal_surf < particle radius for t in [0, dt].
  // Using for x(t) = x_oldPos + t * deltaX/dt and solving for lambda=t/dt
  MFloat lambda = (nx0 - na + radius) / nb;

  // 1) check if a wall-collision happened within the timestep (lambda < 1.0)
  //    meaning that the oldPosition and position lie on different sides of the wall
  //    * lambda > 1.0 means the collision will not happen in this timestep
  //    * lambda < 0.0 means the particle is outside might need to be pushed back inside
  //    * nb > 0.0, means particle is already moving back in
  if(lambda > 1.0 || nb > 0.0) {
    return;
  }

  const MFloat calculatedLambda = lambda;
  // If collision happened already, limit lamba= to -1 (i.e. one timestep)
  if(lambda < -1.0) lambda = -1.0;

  // position of the particle-wall collision
  std::array<MFloat, nDim> oldC{};
  for(MInt i = 0; i < nDim; i++) {
    oldC[i] = m_oldPos[i] + lambda * deltaX[i];
  }

  // considered particle velocity for the wall-collision step
  std::array<MFloat, nDim> usedVel{};
  MFloat velMag = 0;
  for(MInt i = 0; i < nDim; i++) {
    usedVel[i] = deltaX[i] / dt;
    velMag += POW2(usedVel[i]);
  }
  velMag = sqrt(velMag);

#ifdef LPT_DEBUG
  if(m_partId == debugPartId) {
    if(velMag > 0.65) {
      std::cerr << "Large velocity!" << velMag << " " << m_partId << " " << m_position[0] << " " << m_position[1] << " "
                << m_position[2] << " " << m_velocity[0] << " " << m_velocity[1] << " " << m_velocity[2] << " "
                << radius << std::endl;
    }
  }
#endif

  // 2) compute the basis w.r.t. the boundary surface and its inverse
  //    for the coordinate transformations
  std::array<MFloat, nDim> n{};
  std::array<MFloat, nDim> u{};
  std::array<MFloat, nDim> v{};
  for(MInt i = 0; i < nDim; i++) {
    n[i] = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i];
    u[i] = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_planeCoordinates[1][i]
           - s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_planeCoordinates[0][i];
    v[i] = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_planeCoordinates[2][i]
           - s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_planeCoordinates[0][i];
  }
  // normalize the basis vector (normal vector n is already normalized)
  maia::math::normalize(u.data(), nDim);
  maia::math::normalize(v.data(), nDim);
  // gram-schmidt orthogonalization of u and v
  MFloat dotProductUV = std::inner_product(u.begin(), u.end(), v.begin(), 0.0);
  for(MInt i = 0; i < nDim; ++i) {
    v[i] -= dotProductUV * u[i];
  }
  // normalize v again
  maia::math::normalize(v.data(), nDim);

  // construct the basis matrix and its inverse
  std::array<std::array<MFloat, nDim>, nDim> bndryBasis{};
  MFloatScratchSpace inverse(nDim, nDim, AT_, "inverse");
  for(MInt i = 0; i < nDim; i++) {
    bndryBasis[i][0] = n[i];
    bndryBasis[i][1] = u[i];
    bndryBasis[i][2] = v[i];
    for(MInt j = 0; j < nDim; j++) {
      inverse(i, j) = bndryBasis[i][j];
    }
  }
  // NOTE: inverse might not be normalized correctly, i.e. magnitude is not conserved!
  maia::math::invert(&inverse(0, 0), 3, 3);

#ifdef LPT_DEBUG
  if(m_partId == debugPartId) {
    std::cerr << "Basis: " << std::endl;
    std::cerr << "       " << bndryBasis[0][0] << " " << bndryBasis[0][1] << " " << bndryBasis[0][2] << std::endl;
    std::cerr << "       " << bndryBasis[1][0] << " " << bndryBasis[1][1] << " " << bndryBasis[1][2] << std::endl;
    std::cerr << "       " << bndryBasis[2][0] << " " << bndryBasis[2][1] << " " << bndryBasis[2][2] << std::endl;

    std::cerr << "Invert: " << std::endl;
    std::cerr << "       " << inverse(0, 0) << " " << inverse(0, 1) << " " << inverse(0, 2) << std::endl;
    std::cerr << "       " << inverse(1, 0) << " " << inverse(1, 1) << " " << inverse(1, 2) << std::endl;
    std::cerr << "       " << inverse(2, 0) << " " << inverse(2, 1) << " " << inverse(2, 2) << std::endl;

    std::array<MFloat, nDim> a1{};
    std::array<MFloat, nDim> a2{};
    std::array<MFloat, nDim> a3{};

    for(MInt i = 0; i < nDim; i++) {
      a1[i] = inverse(i, 0);
      a2[i] = inverse(i, 1);
      a3[i] = inverse(i, 2);
    }

    maia::math::normalize(&a1[0], nDim);
    maia::math::normalize(&a2[0], nDim);
    maia::math::normalize(&a3[0], nDim);


    MFloatScratchSpace inverse2(nDim, nDim, AT_, "inverse");
    // for(MInt i = 0; i < nDim; i++){
    //  for(MInt j = 0; j < nDim; j++){
    //    inverse2(i,j) = inverse(i,j);
    //  }
    //}

    for(MInt i = 0; i < nDim; i++) {
      inverse2(i, 0) = a1[i];
      inverse2(i, 1) = a2[i];
      inverse2(i, 2) = a3[i];
    }

    std::cerr << "Invert2:" << std::endl;
    std::cerr << "       " << inverse2(0, 0) << " " << inverse2(0, 1) << " " << inverse2(0, 2) << std::endl;
    std::cerr << "       " << inverse2(1, 0) << " " << inverse2(1, 1) << " " << inverse2(1, 2) << std::endl;
    std::cerr << "       " << inverse2(2, 0) << " " << inverse2(2, 1) << " " << inverse2(2, 2) << std::endl;

    maia::math::invert(&inverse2(0, 0), 3, 3);

    std::cerr << "Invert3:" << std::endl;
    std::cerr << "       " << inverse2(0, 0) << " " << inverse2(0, 1) << " " << inverse2(0, 2) << std::endl;
    std::cerr << "       " << inverse2(1, 0) << " " << inverse2(1, 1) << " " << inverse2(1, 2) << std::endl;
    std::cerr << "       " << inverse2(2, 0) << " " << inverse2(2, 1) << " " << inverse2(2, 2) << std::endl;
  }
#endif

  // 3) transform the particle velocity into the boundary coordinate system (forward transformation)
  std::array<MFloat, nDim> localVel{};
  for(MInt i = 0; i < nDim; i++) {
    localVel[i] = 0;
  }
  for(MInt i = 0; i < nDim; i++) {
    for(MInt j = 0; j < nDim; j++) {
      localVel[i] += inverse(i, j) * usedVel[j];
    }
  }

  // 4) Compute rebound velocity by mirrowing the component normal to the boundary
  //    and adding the slip velocity of the wall
  std::array<MFloat, nDim> localReboundVel{};
  localReboundVel[0] = -localVel[0];
  for(MInt i = 1; i < nDim; i++) {
    localReboundVel[i] = localVel[i];
  }

#ifdef LPT_DEBUG
  if(m_partId == debugPartId) {
    std::cerr << "RB " << localReboundVel[0] << " " << localReboundVel[1] << " " << localReboundVel[nDim - 1] << " "
              << velMag << std::endl;
  }
#endif

  // 5) Backwards transformation of the rebound velocity into the x,y,z coordinate system
  std::array<MFloat, nDim> reboundVel{};
  for(MInt i = 0; i < nDim; i++) {
    reboundVel[i] = 0.0;
    for(MInt j = 0; j < nDim; j++) {
      reboundVel[i] += bndryBasis[i][j] * localReboundVel[j];
    }
  }
  maia::math::normalize(reboundVel.data(), nDim);
  for(MInt i = 0; i < nDim; i++) {
    reboundVel[i] *= velMag;
  }
  for(MInt i = 0; i < nDim; i++) {
    m_oldVel[i] = reboundVel[i] + bodyVel[i];
    m_velocity[i] = reboundVel[i] + bodyVel[i];
    m_accel[i] = 0.0;
  }

  // 6) Computation of the new particle position
  if(lambda > 0.0) { // and < 1.0
    for(MInt i = 0; i < nDim; i++) {
      m_position[i] = oldC[i] + (1.0 - lambda) * dt * m_oldVel[i];
    }
  } else {
    // Particle is outside domain or at least overlapping with the surface
    // --> move particle back to the surface using translation in normal direction
    MFloat normalDistance = F0;
    for(MInt i = 0; i < nDim; i++) {
      MFloat surfNormal = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i];
      MFloat surfCoord = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_coordinates[i];
      normalDistance += surfNormal * (m_oldPos[i] - surfCoord);
    }
    for(MInt i = 0; i < nDim; i++) {
      MFloat surfNormal = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i];
      m_position[i] = m_oldPos[i] + (radius - normalDistance) * surfNormal;
    }
  }

  // 7) update the new cellId and check if the particle should be communicated
  this->initProperties();
  hadWallColl() = true;
  checkCellChange(m_oldPos.data());

  // if the push-back is over-extended and the particle moves out-side of an opposide surface/cell
  // the push-back is reduced to the surface
  if(lambda <= 0.0 && (m_cellId < 0 || m_cellId > s_backPtr->a_noCells() - 1 || !s_backPtr->a_isValidCell(m_cellId))) {
    m_cellId = lastValidCellId;
    for(MInt i = 0; i < nDim; i++) {
      m_position[i] = oldC[i];
    }
    this->initProperties();
    hadWallColl() = true;
    checkCellChange(m_oldPos.data());
    std::cerr << "Push back failed for particle " << m_partId << ". Reducing it to surface."
              << " newPos: " << m_cellId << " " << m_position[0] << " " << m_position[1] << " " << m_position[nDim - 1]
              << std::endl;

    // if new position is still invalid but old position was valid, move the particle back to that position
    if(s_backPtr->a_isValidCell(m_oldCellId)
       && (m_cellId < 0 || m_cellId > s_backPtr->a_noCells() - 1 || !s_backPtr->a_isValidCell(m_cellId))) {
      m_cellId = lastValidCellId;
      for(MInt i = 0; i < nDim; i++) {
        m_position[i] = m_oldPos[i];
      }
      this->initProperties();
      hadWallColl() = true;
      checkCellChange(m_oldPos.data());
      std::cerr << "Push back failed again for particle " << m_partId << ". Resetting particle to old position."
                << " newPos: " << m_cellId << " " << m_position[0] << " " << m_position[1] << " "
                << m_position[nDim - 1] << std::endl;

      if(m_cellId < 0 || m_cellId > s_backPtr->a_noCells() - 1 || !s_backPtr->a_isValidCell(m_cellId)) {
        std::cerr << "Push back for particle " << m_partId << " kept failing. Particle could not be rescued."
                  << std::endl;
      }
    }
  }

  // 8) compute old particle coordinate, based on backwards movement from the collision
  //    meaning the position, the particle would have had with current velocity and
  //    without wall interaction
  for(MInt i = 0; i < nDim; i++) {
    if(lambda > 0.0) {
      m_oldPos[i] = oldC[i] + (lambda - 1.0) * dt * m_oldVel[i];
    } else {
      m_oldPos[i] = m_position[i] - dt * m_oldVel[i];
    }
  }

#ifdef LPT_DEBUG
  if(m_partId == debugPartId) {
    std::cerr << "AC " << m_velocity[0] << " " << m_velocity[1] << " " << m_velocity[2] << std::endl;
  }
#endif

  if(m_cellId < 0 || m_cellId > s_backPtr->a_noCells() - 1 || !s_backPtr->a_isValidCell(m_cellId)) {
    std::cerr << "Particle at invalid cell after collision " << m_cellId << std::endl;
    if(m_cellId > 0 && m_cellId < s_backPtr->a_noCells()) {
      std::cerr << s_backPtr->a_isValidCell(m_cellId) << std::endl;
    }
    std::cerr << " " << lastValidCellId << " " << s_backPtr->a_isValidCell(m_oldCellId) << " " << calculatedLambda
              << " newPos " << m_position[0] << " " << m_position[1] << " " << m_position[nDim - 1] << " "
              << s_backPtr->domainId() << " " << m_partId << std::endl;
    std::cerr << "new-OldPos " << m_oldPos[0] << " " << m_oldPos[1] << " " << m_oldPos[nDim - 1] << std::endl;
    std::cerr << "Surface " << s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_coordinates[0] << " "
              << s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_coordinates[1] << " "
              << s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_coordinates[nDim - 1] << std::endl;
    std::cerr << "S-n " << s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[0] << " "
              << s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[1] << " "
              << s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[nDim - 1] << std::endl;
    std::cerr << "CollP " << oldC[0] << " " << oldC[1] << " " << oldC[nDim - 1] << std::endl;
    std::cerr << " deltaX " << deltaX[0] << " " << deltaX[1] << " " << deltaX[nDim - 1] << std::endl;
    std::cerr << "Vel1 " << usedVel[0] << " " << usedVel[1] << " " << usedVel[nDim - 1] << std::endl;
    std::cerr << "Vel2 " << m_oldVel[0] << " " << m_oldVel[1] << " " << m_oldVel[nDim - 1] << std::endl;
    std::cerr << "RbVel " << reboundVel[0] << " " << reboundVel[1] << " " << reboundVel[nDim - 1] << std::endl;
    isInvalid() = true;
    reqSend() = false;
  } else {
    // ensure communication even-though the particle was already in a halo-cell
    if(s_backPtr->a_isHalo(m_cellId)) {
      reqSend() = true;
      isWindow() = false;
    }
  }
}

/// \brief wall-particle collision step
/// \author Tim Wegmann
/// \date   Sep 2022
template <MInt nDim>
void LPTBase<nDim>::wallParticleCollision() {
  const MFloat dt = s_backPtr->m_timeStep;

  const MInt bndId = s_backPtr->a_bndryCellId(m_cellId) > -1 ? s_backPtr->a_bndryCellId(m_cellId)
                                                             : s_backPtr->a_bndryCellId(m_oldCellId);

  if(bndId < 0) return;
  if(isInvalid()) return;
  if(hadWallColl()) return;

  MFloat bodyVelMag = 0;
  MFloat nBodyVel = 0;
  for(MInt i = 0; i < nDim; i++) {
    bodyVelMag += POW2(s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_velocity[i]);
    nBodyVel += s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_velocity[i]
                * s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i];
  }
  bodyVelMag = sqrt(bodyVelMag);

  // stationary or out-ward moving body
  // body-velocity is only added when moving inward i.e. during compression stroke!
  if(bodyVelMag < MFloatEps || nBodyVel < 0) return;

  MFloat nx0 = 0.0;

  std::array<MFloat, nDim> deltaS{};
  std::array<MFloat, nDim> vecSx{};
  for(MInt i = 0; i < nDim; i++) {
    deltaS[i] = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_velocity[i] * dt;
    vecSx[i] = m_position[i] - s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_coordinates[i];
  }

  maia::math::normalize(vecSx.data(), nDim);
  for(MInt i = 0; i < nDim; i++) {
    nx0 += s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_normal[i] * vecSx[i];
  }

  if(nx0 > 0) return;

  // simply add surface velocity to the particle and apply shift
  for(MInt i = 0; i < nDim; i++) {
    m_oldPos[i] = m_position[i];
    m_position[i] = s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_coordinates[i] + deltaS[i];
    m_velocity[i] = m_velocity[i] + s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_velocity[i];
    m_oldVel[i] = m_oldVel[i] + s_backPtr->m_bndryCells->a[bndId].m_srfcs[0]->m_velocity[i];
  }

  // 6) update the new cellId and check if the particle should be communicated
  this->initProperties();
  hadWallColl() = true;
  checkCellChange(m_oldPos.data());

  if(m_cellId < 0 || m_cellId > s_backPtr->a_noCells() - 1 || !s_backPtr->a_isValidCell(m_cellId)) {
    std::cerr << "Particle at invalid cell after collision-2 " << m_cellId << std::endl;
    isInvalid() = true;
    reqSend() = false;
  } else {
    // ensure communication even-though the particle was already in a halo-cell
    if(s_backPtr->a_isHalo(m_cellId)) {
      reqSend() = true;
      isWindow() = false;
    }
  }
}

// Explicit instantiations for 2D and 3D
template class LPTBase<2>;
template class LPTBase<3>;
