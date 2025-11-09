// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HILBERT_H_
#define HILBERT_H_

#include <algorithm>
#include <array>
#include <climits>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"

namespace maia {
namespace grid {
namespace hilbert {

// Detail_ namespace contains implementation details and should never be used
// directly from outside this namespace. For the actual Hilbert index function
// and comments refer to the method "index(...)" below.
namespace detail_ {

// Empty class template to produce error if not specialized below
template <MInt nDim, typename FloatType, typename IdType>
struct Impl;


// Partial specialization for 2-dimensional case (for info on algorithm, see
// "index(...)" method below)
template <typename FloatType, typename IdType>
struct Impl<2, FloatType, IdType> {
  static IdType f(const FloatType* const coordinates, const IdType level) {
    // Define dimension variable for easy reference
    const IdType nDim = 2;

    // Reset index and save coordinates into working array that will be used in
    // the algorithm
    IdType index = 0;
    std::array<FloatType, nDim> x;
    std::copy_n(coordinates, nDim, x.begin());

    // Calculate contribution to Hilbert index iteratively for each level
    for(IdType l = 0; l < level; l++) {
      // SX and SY determine the quadrant in which the coordinates are
      const IdType SX = x[0] < 0.5 ? 0 : 1;
      const IdType SY = x[1] < 0.5 ? 0 : 1;

      // Sn is 1 if coordinates are in quadrant n and 0 otherwise
      // Some quadrants n,m are treated equally and a Snm variable is used
      const IdType S0 = (1 - SX) * (1 - SY);
      const IdType S12 = SY;
      const IdType S3 = SX * (1 - SY);

      // Multiplier is the weight by which a contribution to the Hilbert index
      // is scaled. Higher levels have higher weights and thus a higher impact
      // on the final index.
      const IdType multiplier = IPOW2(nDim * (level - 1 - l));

      // Update index by determining the local quadrand and scaling it
      index += ((1 + SX) * S12 + 3 * S3) * multiplier;

      // Finally, new transformed (mirrored, scaled and translated) coordinates
      // are determined for the next step in the algorithm
      std::array<FloatType, nDim> transformed;
      transformed[0] = 2.0 * x[1] * S0 + (2.0 * x[0] - SX) * S12 + (-2.0 * x[1] + 1.0) * S3;
      transformed[1] = 2.0 * x[0] * S0 + (2.0 * x[1] - 1.0) * S12 + 2.0 * (-x[0] + 1.0) * S3;
      x = transformed;
    }

    return index;
  }
};


// Partial specialization for 3-dimensional case (for info on algorithm, see
// "index(...)" method below)
template <typename FloatType, typename IdType>
struct Impl<3, FloatType, IdType> {
  static IdType f(const FloatType* const coordinates, const IdType level) {
    // Define dimension variable for easy reference
    const IdType nDim = 3;

    // Reset index and save coordinates into working array that will be used in
    // the algorithm
    IdType index = 0;
    std::array<FloatType, nDim> x;
    std::copy_n(coordinates, nDim, x.begin());

    // Calculate contribution to Hilbert index iteratively for each level
    for(IdType l = 0; l < level; l++) {
      // SX and SY determine the quadrant in which the coordinates are
      const IdType SX = x[0] < 0.5 ? 0 : 1;
      const IdType SY = x[1] < 0.5 ? 0 : 1;
      const IdType SZ = x[2] < 0.5 ? 0 : 1;

      // Sn is 1 if coordinates are in quadrant n and 0 otherwise
      // Some quadrants n,m are treated equally and a Snm variable is used
      const IdType S0 = (1 - SX) * (1 - SY) * (1 - SZ);
      const IdType S12 = SY * (1 - SZ);
      const IdType S34 = SX * (1 - SY);
      const IdType S56 = SY * SZ;
      const IdType S7 = (1 - SX) * (1 - SY) * SZ;

      // Multiplier is the weight by which a contribution to the Hilbert index
      // is scaled. Higher levels have higher weights and thus a higher impact
      // on the final index.
      const IdType multiplier = IPOW2(nDim * (level - 1 - l));

      // Update index by determining the local quadrand and scaling it
      index += ((1 + SX) * S12 + (3 + SZ) * S34 + (6 - SX) * S56 + 7 * S7) * multiplier;

      // Finally, new transformed (mirrored, scaled and translated) coordinates
      // are determined for the next step in the algorithm
      std::array<FloatType, nDim> transformed;
      transformed[0] =
          2.0 * x[0] * S0 + 2.0 * x[2] * S12 + (-2.0 * x[1] + 1.0) * S34 + 2.0 * (-x[2] + 1.0) * S56 + 2.0 * x[0] * S7;
      transformed[1] = 2.0 * x[2] * S0 + (2.0 * x[1] - 1.0) * S12 + 2.0 * (-x[0] + 1.0) * S34 + (2.0 * x[1] - 1.0) * S56
                       + 2.0 * (-x[2] + 1.0) * S7;
      transformed[2] = 2.0 * x[1] * S0 + (2.0 * x[0] - SX) * S12 + (2.0 * x[2] - SZ) * S34
                       + (-2.0 * x[0] + 1.0 + SX) * S56 + (-2.0 * x[1] + 1.0) * S7;
      x = transformed;
    }

    return index;
  }
};

} // namespace detail_


/// \brief Return Hilbert index for given location and level in 2D or 3D.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Peter Philippen
/// \date 2014-12-08
///
/// \tparam nDim Dimension for which the Hilbert index should be generated.
/// \tparam FloatType The floating point type for the coordinates.
/// \tparam IdType The integer type that should be used for the index.
/// \param[in] x Pointer to nDim coordinates for which the index is calculated.
///              The coordinates must be within the unit square (2D,
///              [0,1]x[0,1])/unit cube (3D, [0,1]x[0,1]x[0,1]).
/// \param[in] level Refinement level (*not* Hilbert level) for which the index
///                  should be calculated.
///
/// \return The calculated Hilbert index.
///
/// The algorithms used here are taken from
///
///     Peter Philippen. Investigation of two different domain decompositioning
///     methods for a cartesian flow solver. Studienarbeit, 2009.
///
/// The original implementation is by Peter, it was only simplified (no
/// recursion!) and its robustness against bad input parameters increased. Also,
/// the original implementation did not have any comments.
///
/// This is the basic idea:
/// - determine the quadrant/octant of the unit hypercube in which the point is
///   located
/// - calculate its Hilbert index and multiply with the level weight
/// - mirror/scale/translate quadrant/octant to new unit hypercube
/// - repeat from beginning
template <MInt nDim, typename FloatType, typename IdType>
IdType index(const FloatType* const x, const IdType level) {
  // Check that level is at least 0
  ASSERT(level >= 0, "Level must be at least zero");

  // Check for integer overflow ("- 1" to account for signed integers, cast to
  // avoid "unsigned/signed comparison" warning)
  ASSERT(nDim * level <= CHAR_BIT * static_cast<IdType>(sizeof(IdType)) - 1, "Integer overflow: level too large");

  // Assert that coordinates are in unit hypercube
  ASSERT(x[0] >= 0.0 && x[0] <= 1.0, "Coordinate 0 out of bounds.");
  ASSERT(x[1] >= 0.0 && x[1] <= 1.0, "Coordinate 1 out of bounds.");
  IF_CONSTEXPR(nDim == 3) { ASSERT(x[2] >= 0.0 && x[2] <= 1.0, "Coordinate 2 out of bounds."); }

  // Call template function
  return detail_::Impl<nDim, FloatType, IdType>::f(x, level);
}


/// \brief
///
/// \author Lennart
/// \date
template <MInt nDim, typename FloatType, typename IdType>
MBool coordinatesToTreeId(IdType& treeId,
                          const FloatType* const x,
                          const IdType level,
                          const FloatType* const centerOfGravity,
                          FloatType const lengthOnLevel0) {
  // Check that level is at least 0
  ASSERT(level >= 0, "Level must be at least zero");

  // Check for integer overflow ("- 1" to account for signed integers, cast to
  // avoid "unsigned/signed comparison" warning)
  ASSERT(nDim * level <= CHAR_BIT * static_cast<IdType>(sizeof(IdType)) - 1, "Integer overflow: level too large");

  MBool exists = true;
  treeId = 0;
  for(IdType k = 0; k < nDim; k++) {
    if(x[k] < centerOfGravity[k] - 0.5 * lengthOnLevel0 || x[k] > centerOfGravity[k] + 0.5 * lengthOnLevel0) {
      exists = false;
    }
  }
  if(!exists) {
    treeId = std::numeric_limits<IdType>::is_signed ? -1 : ~0;
  } else {
    IdType bitCount = 0;
    IdType coord[3];
    constexpr IdType id2 = (IdType)2;
    for(IdType k = 0; k < nDim; k++) {
      coord[k] = (IdType)((x[k] - centerOfGravity[k] + 0.5 * lengthOnLevel0) * FPOW2(level) / lengthOnLevel0);
    }
    for(IdType lvl = level - 1; lvl >= 0; lvl--) {
      for(IdType k = 0; k < nDim; k++) {
        IdType bit = (coord[k] / (IdType)IPOW2(lvl)) % id2;
        treeId |= bit << bitCount;
        bitCount++;
      }
    }
  }
  return exists;
}

/// \brief
///
/// \author Lennart
/// \date
template <MInt nDim, typename FloatType, typename IdType>
void treeIdToCoordinates(FloatType* const x,
                         const IdType treeId,
                         const IdType level,
                         const FloatType* const centerOfGravity,
                         const FloatType lengthOnLevel0) {
  // Check that level is at least 0
  ASSERT(level >= 0, "Level must be at least zero");

  // Check for integer overflow ("- 1" to account for signed integers, cast to
  // avoid "unsigned/signed comparison" warning)
  ASSERT(nDim * level <= CHAR_BIT * static_cast<IdType>(sizeof(IdType)) - 1, "Integer overflow: level too large");

  IdType bitCount = 0;
  constexpr IdType id1 = (IdType)1;
  for(IdType j = 0; j < nDim; j++) {
    x[j] = centerOfGravity[j];
  }
  for(IdType lvl = 1; lvl <= level; lvl++) {
    for(IdType j = 0; j < nDim; j++) {
      IdType tmpBit = (treeId >> bitCount) & id1;
      x[j] += (tmpBit ? 1.0 : -1.0) * lengthOnLevel0 * FFPOW2(lvl + 1);
      bitCount++;
    }
  }
}

} // namespace hilbert
} // namespace grid
} // namespace maia

#endif // HILBERT_H_
