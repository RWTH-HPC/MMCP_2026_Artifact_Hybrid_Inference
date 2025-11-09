// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGGALERKINPROJECTION_H_
#define DGGALERKINPROJECTION_H_

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "dgcartesianinterpolation.h"
#include "enums.h"
#include "sbpcartesianinterpolation.h"

// Forward declare auxiliary method
namespace maia {
namespace dg {
namespace projection {
template <MInt nDim>
inline void indices(const MInt i, const MInt n, MInt* const ids);
}
} // namespace dg
} // namespace maia

/// \brief Performs the Galerkin projection for a given polynomial degree.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-07-17
template <MInt nDim>
class DgGalerkinProjection {
  // Member methods
 public:
  DgGalerkinProjection() = default;
  DgGalerkinProjection(const MInt polyDeg, const MInt noLevels, const MFloat lengthLevel0, const MInt solverId);
  ~DgGalerkinProjection();

  void calcProjectionInformation(const MInt noDonorCells,
                                 const MInt* donorCellLvls,
                                 const MFloat* donorCellCoordinates,
                                 const MInt targetLvl,
                                 const MFloat* targetElementCenter,
                                 std::vector<MInt>& donorIds,
                                 std::vector<MInt>& donorLvls) const;

  void apply(const MInt noDonorCells, const MInt noNonMappedCells, const MInt* donorCellIds, const MInt* donorCellLvls,
             const MInt* donorCellPos, const MFloat* const donorField, const MFloat* const defaultValues,
             const MInt noVars, MFloat* const targetField) const;

  void calcConservationError(const MInt noDonorCells, const MInt* donorCellIds, const MInt* donorCellLvls,
                             const MFloat* donorField, const MInt noVars, const MFloat targetLength,
                             const MFloat* targetField, MFloat* errors) const;

  void calcL2Error(const MInt noDonorCells, const MInt* donorCellIds, const MInt* donorCellLvls,
                   const MInt* donorCellPos, const MFloat* donorField, const MInt noVars, const MFloat targetLength,
                   const MFloat* targetField, MFloat* errors) const;

 private:
  void init(MInt polyDeg, const MInt noLevels, const MFloat lengthLevel0, const MInt solverId);
  void calcInterpolationMatrices();

  // Member variables
 private:
  // Polynomial degree.
  MInt m_polyDeg = -1;
  // Number of nodes
  MInt m_noNodesXD = -1;
  // Relative number of refinement levels (i.e. range between the coarsest and
  // the finest level)
  MInt m_noLevels = -1;
  // Cell lenght at level zero
  MFloat m_lengthLevel0 = -1.0;

  // Solver id (needed for reading properties)
  MInt m_solverId = -1;

  // Integration weights
  MFloatTensor m_wInt{};

  // Galerkin projection vectors
  // For each level: store the mixed mass matrix entry of each donor cell on
  // this level for each target node
  std::vector<MFloatTensor> m_MTD{};

  // Evaluated Lagrange polynomials for the projection error calculation
  std::vector<MFloatTensor> m_polynomials{};
  // Tolerance for conservation error calculation of the projection
  const MFloat m_tolerance = 1e-12;

  // Cell center positions for all cells on all considered levels
  MFloatTensor m_cellCenters{};
  // Cell lengths for all considered levels
  MFloatTensor m_cellLengths{};
};


/// \brief Apply the Galerkin projection using the provided information.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
///
/// \param[in] noDonorCells Number of donor cells.
/// \param[in] noNonMappedCells Number of non mapped cells/volumes.
/// \param[in] donorCellIds Ids of the donor cells involved in the projection.
/// \param[in] donorCellLevels Relative levels of the donor cells (including non-mapped cells).
/// \param[in] donorCellPos Donor cell positions (index on cell level; incl. non-mapped cells).
/// \param[in] donorField Variable values of all donor cells.
/// \param[in] noVars Number of variables.
/// \param[out] targetField Target element variables.
template <MInt nDim>
void DgGalerkinProjection<nDim>::apply(const MInt noDonorCells,
                                       const MInt noNonMappedCells,
                                       const MInt* donorCellIds,
                                       const MInt* donorCellLevels,
                                       const MInt* donorCellPos,
                                       const MFloat* const donorField,
                                       const MFloat* const defaultValues,
                                       const MInt noVars,
                                       MFloat* const targetField) const {
  MFloatTensor target(targetField, m_noNodesXD, noVars);
  // Set target field to zero
  std::fill_n(&targetField[0], m_noNodesXD * noVars, 0.0);

  // Compute target field
  // Loop over all donor cells
  for(MInt l = 0; l < noDonorCells; l++) {
    // Weights of this donor cell on all target nodes
    // TODO labels:DG,totest,toenhance Check at some point in the future if the very indirect access to
    //       the projection weights is too slow and should be optimized
    const MFloat* const weight = &m_MTD[donorCellLevels[l]](donorCellPos[l], 0);
    // Variables of donor cell
    const MFloat* const donorVars = &donorField[donorCellIds[l] * noVars];

    // Loop over all target nodes
    for(MInt i = 0; i < m_noNodesXD; i++) {
      // Loop over all variables
      for(MInt v = 0; v < noVars; v++) {
        target(i, v) += weight[i] * donorVars[v];
      }
    }
  }

  // Use defaul values for all non-mapped cells (i.e. cells of the mapped volume that are missing)
  for(MInt l = noDonorCells; l < noDonorCells + noNonMappedCells; l++) {
    const MFloat* const weight = &m_MTD[donorCellLevels[l]](donorCellPos[l], 0);
    // Loop over all target nodes
    for(MInt i = 0; i < m_noNodesXD; i++) {
      // Loop over all variables
      for(MInt v = 0; v < noVars; v++) {
        target(i, v) += weight[i] * defaultValues[v];
      }
    }
  }
}


/// \brief Calculate the Galerkin projection conservation error and check
///        against tolerance.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
///
/// \param[in] noDonorCells Number of donor cells.
/// \param[in] donorCellLevels Relative levels of donor cells.
/// \param[in] donorField Donor field values.
/// \param[in] noVars Number of variables.
/// \param[in] targetLength Length of the target element.
/// \param[in] targetField Computed target field.
/// \param[out] errors Computed conservation errors.
///
/// The conservation error of the Galerkin projection is defined as:
/// E_C = |\int_{Q_T} (q_D - q_T) dV|
/// i.e. the absolute difference of the donor field and the target field
/// integrated over the target element.
/// If any of the computed conservation errors exceeds a tolerance of 1e-12 the
/// simulation will abort.
template <MInt nDim>
void DgGalerkinProjection<nDim>::calcConservationError(const MInt noDonorCells, const MInt* donorCellIds,
                                                       const MInt* donorCellLevels, const MFloat* donorField,
                                                       const MInt noVars, const MFloat targetLength,
                                                       const MFloat* targetField, MFloat* errors) const {
  TRACE();

  const MInt noNodes1D = m_polyDeg + 1;
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;

  const MFloatTensor target(const_cast<MFloat*>(targetField), noNodes1D, noNodes1D, noNodes1D3, noVars);

  // Jacobian of the element transformation
  const MFloat elementJacobian = (nDim == 2) ? POW2(F1B2 * targetLength) : POW3(F1B2 * targetLength);

  MFloatScratchSpace donorIntegral(noVars, AT_, "donorIntegral");
  std::fill(donorIntegral.begin(), donorIntegral.end(), 0.0);
  MFloatScratchSpace targetIntegral(noVars, AT_, "targetIntegral");
  std::fill(targetIntegral.begin(), targetIntegral.end(), 0.0);

  // Calculate donor integral over all donor cells ...
  for(MInt l = 0; l < noDonorCells; l++) {
    const MInt dLevel = donorCellLevels[l];
    const MInt nCellsOnLevel = IPOW2(dLevel);

    const MFloat donorLength = targetLength / nCellsOnLevel;
    const MFloat volume = (nDim == 2) ? POW2(donorLength) : POW3(donorLength);

    const MInt donorCellOffset = donorCellIds[l] * noVars;
    // ... for all variables
    for(MInt varId = 0; varId < noVars; varId++) {
      donorIntegral[varId] += volume * donorField[donorCellOffset + varId];
    }
  }

  // Calculate target integral for all variables
  IF_CONSTEXPR(nDim == 2) { // 2D
    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D; j++) {
        const MFloat weight = m_wInt(i) * m_wInt(j);
        for(MInt varId = 0; varId < noVars; varId++) {
          targetIntegral[varId] += weight * target(i, j, 0, varId);
        }
      }
    }
  }
  else { // 3D
    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt k = 0; k < noNodes1D; k++) {
          const MFloat weight = m_wInt(i) * m_wInt(j) * m_wInt(k);
          for(MInt varId = 0; varId < noVars; varId++) {
            targetIntegral[varId] += weight * target(i, j, k, varId);
          }
        }
      }
    }
  }

  // Multiply target integral with jacobian of the transformation
  for(MInt varId = 0; varId < noVars; varId++) {
    targetIntegral[varId] *= elementJacobian;
  }

  // Check if donor and target integrals match
  for(MInt varId = 0; varId < noVars; varId++) {
    errors[varId] = ABS(donorIntegral[varId] - targetIntegral[varId]);

    // Exit with error if conservation error is too large
    if(errors[varId] > m_tolerance) {
      std::stringstream errorMsg;
      errorMsg << "targetIntegral doesn't match donorIntegral! (donorIntegral = " << std::scientific
               << donorIntegral[varId] << "; targetIntegral = " << targetIntegral[varId]
               << "; error = " << errors[varId] << "; m_tolerance = " << m_tolerance << ")";

      TERMM(1, errorMsg.str());
    }
  }
}


/// \brief Calculate the L2 projection error.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
///
/// \param[in] noDonorCells Number of donor cells.
/// \param[in] donorCellLevels Relative levels of donor cells.
/// \param[in] donorCellPos Positions of donor cells.
/// \param[in] donorField Donor field values.
/// \param[in] noVars Number of variables.
/// \param[in] targetLength Length of the target element.
/// \param[in] targetField Computed target field.
/// \param[out] errors Computed projection errors.
///
/// The L2 projection error is defined as:
/// E_L2 = || q_D - q_T ||_2^2 = \int_{Q_T} (q_D - q_T)^2 dV
/// with q_D and q_T the donor and target field. The target field is expanded in
/// terms of its basis and the integration over the target element is written as
/// the integration over all donor cells, since the donor field is piecewise
/// constant.
template <MInt nDim>
void DgGalerkinProjection<nDim>::calcL2Error(const MInt noDonorCells, const MInt* donorCellIds,
                                             const MInt* donorCellLevels, const MInt* donorCellPos,
                                             const MFloat* donorField, const MInt noVars, const MFloat targetLength,
                                             const MFloat* targetField, MFloat* errors) const {
  TRACE();

  using namespace maia::dg::projection;

  const MInt noNodes1D = m_polyDeg + 1;
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;

  const MFloatTensor target(const_cast<MFloat*>(targetField), m_noNodesXD, noVars);
  const MFloat elementJacobian = (nDim == 2) ? POW2(F1B2 * targetLength) : POW3(F1B2 * targetLength);

  MInt intervals[3], index[3];

  std::fill_n(&errors[0], noVars, 0.0);

  // Loop over all donor cells
  for(MInt l = 0; l < noDonorCells; l++) {
    const MInt dLevel = donorCellLevels[l];
    const MInt nCellsOnLevel = IPOW2(dLevel);

    const MFloat donorXiLength = m_cellLengths[dLevel];
    const MFloat donorCellJacobian = (nDim == 2) ? POW2(F1B2 * donorXiLength) : POW3(F1B2 * donorXiLength);

    const MInt donorCellOffset = donorCellIds[l] * noVars;

    // Compute donor cell index for all directions
    indices<nDim>(donorCellPos[l], nCellsOnLevel, intervals);

    // Loop over all integration nodes on donor cell
    for(MInt u = 0; u < noNodes1D; u++) {
      for(MInt k = 0; k < noNodes1D; k++) {
        for(MInt j = 0; j < noNodes1D3; j++) {
          std::vector<MFloat> sum(noVars, 0.0);

          // Loop over all nodes of element
          for(MInt i = 0; i < m_noNodesXD; i++) {
            indices<nDim>(i, noNodes1D, index); // Node position

            // Multiply Lagrange polynomials
            MFloat tmp =
                m_polynomials[dLevel](intervals[0], u, index[0]) * m_polynomials[dLevel](intervals[1], k, index[1]);
            IF_CONSTEXPR(nDim == 3) { tmp *= m_polynomials[dLevel](intervals[2], j, index[2]); }

            // Add contribution of current element node
            for(MInt varId = 0; varId < noVars; varId++) {
              sum[varId] += tmp * target(i, varId);
            }
          }

          MFloat factor = m_wInt(u) * m_wInt(k) * elementJacobian * donorCellJacobian;
          IF_CONSTEXPR(nDim == 3) { factor *= m_wInt(j); }

          // Add to projection error
          for(MInt varId = 0; varId < noVars; varId++) {
            const MFloat square = POW2(donorField[donorCellOffset + varId] - sum[varId]);
            errors[varId] += factor * square;
          }
        }
      }
    }
  }
}


namespace maia {
namespace dg {
namespace projection {

/// \brief Calculate the 2D/3D indices for a given scalar id for accessing a
///        field of n^nDim points.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-06-01
///
/// \param[in] i Scalar index.
/// \param[in] n Number of points in one dimension.
/// \param[out] ids Computed indices for all spatial directions.
///
/// This function is the inverse of, e.g. computing a (square) matrix entry
/// position in memory as i*n+j. The last dimension always varies the fastest.
template <MInt nDim>
inline void indices(const MInt i, const MInt n, MInt* const ids) {
  const MInt f1 = i / n;

  IF_CONSTEXPR(nDim == 2) {
    ids[0] = f1;
    ids[1] = i - f1 * n;
  }
  else {
    const MInt f2 = i / (n * n);

    ids[0] = f2;
    ids[1] = f1 - f2 * n;
    ids[2] = i - f1 * n;
  }
}

} // namespace projection
} // namespace dg
} // namespace maia


#endif /* DGGALERKINPROJECTION_H_ */
