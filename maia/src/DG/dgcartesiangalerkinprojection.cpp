// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "dgcartesiangalerkinprojection.h"

#include <array>
#include "IO/context.h"
#include "MEMORY/scratch.h"
#include "UTIL/tensor.h"
#include "UTIL/timer.h"
#include "property.h"

using namespace maia::dg::interpolation;
using namespace maia::dg::projection;
using namespace std;


/// \brief Constructor passes arguments to init().
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
///
/// \param[in] polyDeg Polynomial degree of projection.
/// \param[in] noLevels Relative number of refinement levels, i.e. the
///                     range between the coarsest and the finest level.
/// \param[in] lenghtLevel0 Cell length on level zero.
/// \param[in] solverId Solver id.
template <MInt nDim>
DgGalerkinProjection<nDim>::DgGalerkinProjection(const MInt polyDeg, const MInt noLevels, const MFloat lengthLevel0,
                                                 const MInt solverId) {
  TRACE();

  init(polyDeg, noLevels, lengthLevel0, solverId);
}


/// \brief Destructor clears all member variables.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
template <MInt nDim>
DgGalerkinProjection<nDim>::~DgGalerkinProjection() {
  TRACE();

  m_polyDeg = -1;
  m_noNodesXD = -1;
  m_noLevels = -1;
  m_lengthLevel0 = -1.0;
  m_solverId = -1;
  m_wInt.clear(), m_MTD.clear();
  m_polynomials.clear();
  m_cellCenters.clear();
  m_cellLengths.clear();
}


/// \brief Initialize the member variables and calculate the interpolation
///        matrices.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
///
/// \param[in] polyDeg Polynomial degree.
/// \parma[in] noLevels Maximum difference in refinement levels between donor
///            and target mesh + 1.
/// \param[in] lenghtLevel0 Cell length on level zero.
/// \param[in] solverId Solver id.
template <MInt nDim>
void DgGalerkinProjection<nDim>::init(const MInt polyDeg,
                                      const MInt noLevels,
                                      const MFloat lengthLevel0,
                                      const MInt solverId) {
  TRACE();

  // Set member variables
  m_polyDeg = polyDeg;
  m_noLevels = noLevels;
  m_lengthLevel0 = lengthLevel0;
  m_noNodesXD = ipow(polyDeg + 1, nDim);
  m_solverId = solverId;

  const MInt noNodes1D = polyDeg + 1;
  const MInt maxNoCellsOnLevel = IPOW2(noLevels - 1);

  // Allocate space for member variables
  m_MTD.resize(m_noLevels);
  m_polynomials.resize(m_noLevels);
  for(MInt level = 0; level < m_noLevels; level++) {
    const MInt noCells1D = IPOW2(level);
    const MInt noCellsOnLevel = ipow(noCells1D, nDim);

    m_MTD[level].resize(noCellsOnLevel, m_noNodesXD);
    m_polynomials[level].resize(noCells1D, noNodes1D, noNodes1D);
  }

  m_cellCenters.resize(noLevels, maxNoCellsOnLevel);
  m_cellLengths.resize(noLevels);

  // Calculate cell center positions and cell lengths for all levels
  for(MInt l = 0; l < noLevels; l++) {
    // Store cell length for this level
    const MInt noCellsOnLevel = IPOW2(l);
    const MFloat cellLength = 2.0 / noCellsOnLevel;
    m_cellLengths[l] = cellLength;

    // Store cell center positions
    for(MInt c = 0; c < noCellsOnLevel; c++) {
      m_cellCenters(l, c) = -1.0 + 0.5 * cellLength + c * cellLength;
    }
  }

  // Calculate the needed interpolation matrices
  calcInterpolationMatrices();
}


/// \brief Calculate and store the needed interpolation matrices.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
///
/// The Galerkin projection can generally be written as:
/// q_T = (M_T)^-1 * sum_{l=1}^{S} (M_TDl * q_Dl)
/// with q_T the target field and q_Dl the piecewise constant solutions
/// of the donor cells 1 to l. The mass matrix M_T contains the integration
/// weights of each node on the diagonal, thus inverting it is trivial.
/// The vector M_TDl contains the weighting of the l-th donor cell solution on
/// each target element node.
/// Given the Cartesian grid structure the vectors M_TDl can be precomputed for
/// all relative refinement levels between the donor and the target grid and all
/// cells on these levels.
/// This method computes the vectors M_TD and premultiplies them with the
/// inverse mass matrix. The projection then has the form:
/// q_T[nodeId] = sum_{l=1}^{S} (MTD[donorLevel][donorPos][nodeId] * q_Dl)
/// where donorLevel and donorPos is the relative level and position of the
/// donor cell relative to the target cell.
///
/// Further, the polynomials are evaluated at all integration nodes on all cells
/// on all levels (in 1D) to allow evaluation of the L2 projection error.
///
/// For more information refer to:
/// "Conservative Galerkin projection for hybrid computational aeroacoustics on
/// hierarchical Cartesian grids", Ansgar Marwege, Bachelor Thesis
template <MInt nDim>
void DgGalerkinProjection<nDim>::calcInterpolationMatrices() {
  TRACE();

  // Check for SBP Mode
  const MBool defaultSbpMode = false;
  const MBool sbpMode = Context::getSolverProperty<MBool>("sbpMode", m_solverId, AT_, &defaultSbpMode);

  // Determine SBP Operator
  const MString defaultSbpOperator = "";
  const MString sbpOperator = Context::getSolverProperty<MString>("sbpOperator", m_solverId, AT_, &defaultSbpOperator);

  MInt initNoNodes1D = -1;
  initNoNodes1D = Context::getSolverProperty<MInt>("initNoNodes", m_solverId, AT_, &initNoNodes1D);

  const MInt noNodes1D = sbpMode ? initNoNodes1D : (m_polyDeg + 1);
  const MInt maxNoCellsOnLevel = IPOW2(m_noLevels - 1);

  // Determine DG integration method (same as in dgsolver.cpp)
  MString defaultDgIntegrationMethod = "DG_INTEGRATE_GAUSS";
  const MInt dgIntegrationMethod = string2enum(
      Context::getSolverProperty<MString>("dgIntegrationMethod", m_solverId, AT_, &defaultDgIntegrationMethod));

  // Determine DG polynomial type (same as in dgsolver.cpp)
  MString defaultDgPolynomialType = "DG_POLY_LEGENDRE";
  const MInt dgPolynomialType =
      string2enum(Context::getSolverProperty<MString>("dgPolynomialType", m_solverId, AT_, &defaultDgPolynomialType));

  // Convert integers to enums
  DgPolynomialType polyType = static_cast<DgPolynomialType>(dgPolynomialType);
  DgIntegrationMethod intMethod = static_cast<DgIntegrationMethod>(dgIntegrationMethod);

  // Create interpolation object
  DgInterpolation interpolation(m_polyDeg, polyType, noNodes1D, intMethod, sbpMode, sbpOperator);

  // Create convenience pointers
  const MFloatVector& wInt = interpolation.m_wInt;
  const MFloatVector& wBary = interpolation.m_wBary;
  const MFloatVector& nodes = interpolation.m_nodes;

  // Store integration weights in class for later use
  m_wInt = wInt;

  // Storage for evaluating the lagrange polynomials at a given position
  ScratchSpace<MFloat> polynomials(m_polyDeg + 1, AT_, "polynomials");
  // Storage for inverse mass matrix (diagonal)
  ScratchSpace<MFloat> MT_inverse(m_noNodesXD, AT_, "MT_inverse");
  // Storage for 1D-vectors
  ScratchSpace<MFloat> vector1D(m_noLevels, maxNoCellsOnLevel, noNodes1D, AT_, "vector1D");

  std::array<MInt, nDim> index;
  // Calculate inverted mass matrix with the integration weights
  IF_CONSTEXPR(nDim == 2) {
    for(MInt i = 0; i < m_noNodesXD; i++) {
      indices<nDim>(i, noNodes1D, &index[0]);
      MT_inverse[i] = 1.0 / (wInt(index[0]) * wInt(index[1]));
    }
  }
  else {
    for(MInt i = 0; i < m_noNodesXD; i++) {
      indices<nDim>(i, noNodes1D, &index[0]);
      MT_inverse[i] = 1.0 / (wInt(index[0]) * wInt(index[1]) * wInt(index[2]));
    }
  }

  // Compute 1D-vectors
  // i.e. on all relative refinement levels: the weight of each donor cell on
  // this level for each target node of the element (in one dimension)
  // Loop over all (relative) refinement levels, start with highest level
  for(MInt level = m_noLevels - 1; level >= 0; level--) {
    const MInt noCellsOnLevel = IPOW2(level);
    const MFloat lengthOnLevel = m_cellLengths[level];

    // Loop over all possible donor cells on this level
    for(MInt l = 0; l < noCellsOnLevel; l++) {
      // Loop over all target nodes
      for(MInt j = 0; j < noNodes1D; j++) {
        MFloat sumXi = 0.0;
        // On the finest level do the integration by quadrature
        if(level == m_noLevels - 1) {
          // Loop over quadrature points
          for(MInt k = 0; k < noNodes1D; k++) {
            // Node position on reference element [-1,1]
            const MFloat xiPos = m_cellCenters(level, l) + 0.5 * lengthOnLevel * nodes[k];

            // Evaluate Lagrange polynomials
            if(sbpMode) {
              calcLinearInterpolationBase(xiPos, noNodes1D, &nodes[0], &polynomials[0]);
            } else {
              calcLagrangeInterpolatingPolynomials(xiPos, m_polyDeg, &nodes[0], &wBary[0], &polynomials[0]);
            }
            // Calculate weight and add to integral
            const MFloat weight = 0.5 * lengthOnLevel * wInt[k];
            sumXi += weight * polynomials[j];
          }
        } else {
          // Combine the integrals of the finer level
          sumXi = vector1D(level + 1, 2 * l, j) + vector1D(level + 1, 2 * l + 1, j);
        }

        // Store projection information
        vector1D(level, l, j) = sumXi;
      }
    }
  }

  // TODO labels:DG,toenhance use symmetry to save on memory and computations?
  MInt cellIndex[MAX_SPACE_DIMENSIONS];
  // Calculate projection vectors M_TD for all cells on all levels (Cartesian
  // product of 1D-vectors)
  for(MInt level = 0; level < m_noLevels; level++) {
    const MInt noCells1D = IPOW2(level);
    const MInt noCellsOnLevel = ipow(noCells1D, nDim);

    // Loop over all target nodes
    for(MInt nodeId = 0; nodeId < m_noNodesXD; nodeId++) {
      // Compute node index in all dimensons from nodeId
      indices<nDim>(nodeId, noNodes1D, &index[0]);

      // Loop over all cells on this level
      for(MInt i = 0; i < noCellsOnLevel; i++) {
        // Compute cell index in all dimensions
        indices<nDim>(i, noCells1D, cellIndex);

        // Assemble projection vector M_TD
        m_MTD[level](i, nodeId) = 1.0;
        for(MInt d = 0; d < nDim; d++) {
          m_MTD[level](i, nodeId) *= vector1D(level, cellIndex[d], index[d]);
        }
        // Premultiply with inverted mass matrix
        m_MTD[level](i, nodeId) *= MT_inverse[nodeId];
      }
    }
  }

  // Evaluate polynomials for error calculation
  for(MInt level = 0; level < m_noLevels; level++) {
    const MInt noCellsOnLevel = IPOW2(level);
    const MFloat lengthOnLevel = m_cellLengths[level];

    // Loop over all possible donor cells on this level
    for(MInt l = 0; l < noCellsOnLevel; l++) {
      // Loop over quadrature points
      for(MInt k = 0; k < noNodes1D; k++) {
        // Node position
        const MFloat xiPos = m_cellCenters(level, l) + 0.5 * lengthOnLevel * nodes[k];

        // Evaluate Lagrange polynomials
        calcLagrangeInterpolatingPolynomials(xiPos, m_polyDeg, &nodes[0], &wBary[0], &m_polynomials[level](l, k, 0));
      }
    }
  }
}


/// \brief Compute the projection information for one target cell.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-07-16
///
/// \param[in] noDonorCells The number of donor cells.
/// \param[in] donorCellLevels The levels of the donor cells.
/// \param[in] donorCellCoordinates Coordinates of the donor cells.
/// \param[in] targetLevel Level of the target cell.
/// \param[in] targetCellCenter Coordinates of the target cell.
/// \param[out] donorIds List of computed donor cell indices.
/// \param[out] donorLevels List of donor cell levels relative to the target
///             cell.
///
/// The projection information for one target cell consists of a list of donor
/// cell indices and a list of donor cell levels. The donor cell level stores
/// the relative level between the donor cell and the target cell. With the
/// donor cell index the position of the donor cell on its level relative to the
/// cell is described, i.e. it is a scalar index instead of nDim-indices for all
/// dimensions.
/// Additionally, the non-mapped volumes are identified and added to the donorIds and donorLevels
/// lists. This means that cells that are missing in the FV-solver since they are inside some
/// geometric object are identified such that for these volumes default values for mapping data to a
/// DG element can be used (without implicitly assuming them all to be zero when these volumes are
/// not handled).
template <MInt nDim>
void DgGalerkinProjection<nDim>::calcProjectionInformation(const MInt noDonorCells, const MInt* donorCellLevels,
                                                           const MFloat* donorCellCoordinates, const MInt targetLevel,
                                                           const MFloat* targetCellCenter, std::vector<MInt>& donorIds,
                                                           std::vector<MInt>& donorLevels) const {
  TRACE();

  // Calculates the position on the unit element corresponding to the cell with
  // given center and width.
  auto positionOnUnitElement = [](MFloat pos, MFloat center, MFloat width) { return 2 * (pos - center) / width; };

  // Define epsilon according to CartesianGrid constructor
  const MFloat eps = 1.0 / FPOW2(30) * m_lengthLevel0;

  const MFloatTensor donorCellCoords(const_cast<MFloat*>(donorCellCoordinates), noDonorCells, nDim);
  const MFloat cellLength = m_lengthLevel0 * FFPOW2(targetLevel);
  MInt intervals[MAX_SPACE_DIMENSIONS];

  IF_CONSTEXPR(nDim == 2) {
    intervals[2] = 0; // Not set in 2D, prevent uninitialized value valgrind errors
  }

  // Determine maximum level difference between donor cells and target cell
  const MInt maxDonorLevel =
      *(std::max_element(&donorCellLevels[0], &donorCellLevels[0] + noDonorCells)) - targetLevel + 1;

  // Initialize storage for mapped volume, i.e. consider each possible cell on all levels
  std::vector<std::vector<MInt>> mappedVolume(maxDonorLevel);
  for(MInt i = 0; i < maxDonorLevel; i++) {
    const MInt noCellsOnLvl = ipow(IPOW2(i), nDim);
    mappedVolume[i].resize(noCellsOnLvl);
    std::fill(mappedVolume[i].begin(), mappedVolume[i].end(), 0);
  }

  // Loop over all donor cells
  for(MInt l = 0; l < noDonorCells; l++) {
    // Relative level of donor cell
    const MInt donorLevel = donorCellLevels[l] - targetLevel;

    const MInt noCells1D = IPOW2(donorLevel);
    const MInt noCells1D3 = (nDim == 3) ? noCells1D : 1;

    // Determine position of donor cell relative to target cell
    for(MInt d = 0; d < nDim; d++) {
      const MFloat posOnUnitElem = positionOnUnitElement(donorCellCoords(l, d), targetCellCenter[d], cellLength);

      intervals[d] = -1;
      // Find interval for current coordinate direction
      for(MInt i = 0; i < noCells1D; i++) {
        const MFloat lowerBound = m_cellCenters(donorLevel, i) - 0.5 * m_cellLengths[donorLevel];
        const MFloat upperBound = m_cellCenters(donorLevel, i) + 0.5 * m_cellLengths[donorLevel];

        // Check if the donor cell center is contained in the current cell
        // interval and if it is 'roughly' the same as the current cell center
        if(lowerBound < posOnUnitElem && posOnUnitElem < upperBound
           && approx(posOnUnitElem, m_cellCenters(donorLevel, i), eps)) {
          intervals[d] = i;
          break;
        }
      }

      if(intervals[d] == -1) {
        stringstream errorMsg;
        errorMsg << "Donor cell interval not found! (level = " << donorLevel << "; position = " << std::scientific
                 << posOnUnitElem << ")";
        TERMM(1, errorMsg.str());
      }
    }

    // Compute scalar donor cell index
    const MInt donorIndex =
        intervals[0] * noCells1D * noCells1D3 + intervals[1] * noCells1D3 + intervals[2] * (nDim == 3);

    // Store donor cell information
    donorIds[l] = donorIndex;
    donorLevels[l] = donorLevel;

    // Keep track of mapped volume, to identify non-mapped/missing cells
    mappedVolume[donorLevel][donorIndex] = 1;

    // Loop down to highest level and mark all descendent cells as mapped
    for(MInt lvl = donorLevel + 1; lvl < maxDonorLevel; lvl++) {
      const MInt noCellsRel1D = IPOW2(lvl - donorLevel);
      const MInt noCellsRel1D3 = (nDim == 3) ? noCellsRel1D : 1;
      const MInt sizeFactor = IPOW2(lvl - donorLevel);

      const MInt noCells1D_ = IPOW2(lvl);
      const MInt noCells1D3_ = (nDim == 3) ? noCells1D_ : 1;

      // Loop over all child positions
      for(MInt i = 0; i < noCellsRel1D; i++) {
        for(MInt j = 0; j < noCellsRel1D; j++) {
          for(MInt k = 0; k < noCellsRel1D3; k++) {
            // Determine child 1D index on its level and mark as mapped
            const MInt volumeIndex = (intervals[0] * sizeFactor + i) * noCells1D_ * noCells1D3_
                                     + (intervals[1] * sizeFactor + j) * noCells1D3_
                                     + (intervals[2] * sizeFactor + k) * (nDim == 3);
            mappedVolume[lvl][volumeIndex] = 1;
          }
        }
      }
    }

    // Loop up to relative level zero and mark coarser cells as mapped
    for(MInt lvl = donorLevel - 1; lvl >= 0; lvl--) {
      const MInt noCells1D_ = IPOW2(lvl);
      const MInt noCells1D3_ = (nDim == 3) ? noCells1D_ : 1;
      const MInt sizeFactor = IPOW2(donorLevel - lvl);

      // Determine 1D index of parent on coarser level and mark as (partially) mapped
      const MInt volumeIndex = std::floor(intervals[0] / sizeFactor) * noCells1D_ * noCells1D3_
                               + std::floor(intervals[1] / sizeFactor) * noCells1D3_
                               + std::floor(intervals[2] / sizeFactor) * (nDim == 3);
      mappedVolume[lvl][volumeIndex] = 1;
    }
  }

  // Append non-mapped cells/volumes to projection information
  // Loop over all donor cell levels from coarse to fine
  for(MInt lvl = 0; lvl < maxDonorLevel; lvl++) {
    const MInt noCellsLvl1D = IPOW2(lvl);
    const MInt noCellsLvl1D3 = (nDim == 3) ? noCellsLvl1D : 1;

    // Loop over all cells on this level
    for(MInt i = 0; i < noCellsLvl1D; i++) {
      for(MInt j = 0; j < noCellsLvl1D; j++) {
        for(MInt k = 0; k < noCellsLvl1D3; k++) {
          // Index of current cell on this level
          const MInt donorIndex = i * noCellsLvl1D * noCellsLvl1D3 + j * noCellsLvl1D3 + k * (nDim == 3);

          // If the current cell on this level is not marked as mapped ...
          if(!mappedVolume[lvl][donorIndex]) {
            // ... append its position and level to the projection information
            donorIds.push_back(donorIndex);
            donorLevels.push_back(lvl);

            std::cerr << globalDomainId() << " non-mapped " << targetCellCenter[0] << " " << lvl << " " << donorIndex
                      << std::endl;
            TERMM(1, "TODO check if this is working correctly");

            // Loop down to highest level and mark all children of this cell as mapped since the
            // full volume of the coarse cell is already mapped
            for(MInt childLvl = lvl + 1; childLvl < maxDonorLevel; childLvl++) {
              const MInt noCellsRel1D = IPOW2(childLvl - lvl);
              const MInt noCellsRel1D3 = (nDim == 3) ? noCellsRel1D : 1;

              // Loop over all cells on this child level and mark as mapped
              for(MInt i2 = 0; i2 < noCellsRel1D; i2++) {
                for(MInt j2 = 0; j2 < noCellsRel1D; j2++) {
                  for(MInt k2 = 0; k2 < noCellsRel1D3; k2++) {
                    const MInt volumeIndex =
                        (i2 + i) * noCellsRel1D * noCellsRel1D3 + (j2 + j) * noCellsRel1D3 + (k2 + k) * (nDim == 3);
                    mappedVolume[childLvl][volumeIndex] = 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Compute total volume of mapped cells and check (should be 1.0)
  MFloat sumVolumes = 0.0;
  for(MUint i = 0; i < donorLevels.size(); i++) {
    const MFloat volume = 1.0 / (ipow(IPOW2(donorLevels[i]), nDim));
    sumVolumes += volume;
  }
  if(!approx(sumVolumes, 1.0, eps)) {
    stringstream errMsg;
    errMsg << "calcProjectionInformation: mapped volume " << std::to_string(sumVolumes) << "(target cell coordinates:";
    for(MInt d = 0; d < nDim; d++) {
      errMsg << " " << targetCellCenter[d];
    }
    errMsg << ")";
    TERMM(1, errMsg.str());
  }
}


template class DgGalerkinProjection<2>;
template class DgGalerkinProjection<3>;
