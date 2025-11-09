// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_PARTICLELIB_H
#define MAIA_PARTICLELIB_H

#include <random>
#include "UTIL/maiamath.h"
#include "globals.h"

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

template <MInt nDim>
class LPTEllipsoidal;

template <MInt nDim>
class LPTSpherical;

template <MInt nDim>
class LPT;

namespace maia {
namespace lpt {

template <MInt nDim>
using ellipsListIterator = typename std::vector<LPTEllipsoidal<nDim>>::iterator;
template <MInt nDim>
using ellipsListIteratorConst = typename std::vector<LPTEllipsoidal<nDim>>::const_iterator;

template <MInt nDim>
using partListIterator = typename std::vector<LPTSpherical<nDim>>::iterator;
template <MInt nDim>
using partListIteratorConst = typename std::vector<LPTSpherical<nDim>>::const_iterator;

template <MInt nDim>
struct subDomainCollector {
  typename std::vector<typename std::vector<LPTSpherical<nDim>>::iterator> subDomain;
};

template <MInt nDim>
struct subDomainCollectorEllipsoid {
  typename std::vector<typename std::vector<LPTEllipsoidal<nDim>>::iterator> subDomain;
};

template <MInt nDim>
struct sendQueueType {
  MLong partPos;
  MInt toCellId;
};

struct partType {
  MInt partId;
  MInt cellId;
  MFloat diam;
  MFloat densRatio;
};

struct partTypeEllipsoid {
  MInt partId;
  MInt cellId;
  MFloat semiMinorAxis;
  MFloat densRatio;
  MFloat aspectRatio;
};


template <class T>
class sort_particleAfterPartIds {
 public:
  // NOTE: with secondary-break-up particle-ids might not be
  // fully unique after a restart anymore, in this case they are additionally
  // sorted by their x-position
  static MBool compare(const T& i, const T& j) {
    if(i.m_partId == j.m_partId) {
      return i.m_position[0] < j.m_position[0];
    }
    return (i.m_partId < j.m_partId);
  }
};

template <class T>
class sort_particleAfterDiameter {
 public:
  static MBool compare(const T& i, const T& j) { return (i.m_diameter < j.m_diameter); }
};

template <class T>
class sortDesc_particleAfterDiameter {
 public:
  static MBool compare(const T& i, const T& j) { return (i.m_diameter > j.m_diameter); }
};

template <class T>
class sort_particleAfterTemperature {
 public:
  static MBool compare(const T& i, const T& j) { return (i.m_temperature > j.m_temperature); }
};

template <class T>
class findPartId {
 public:
  static MBool compare(const T& i, const T& j) { return (i.m_partId < j.m_partId); }
  static MBool compare(const MInt i, const T& j) { return (i < j.m_partId); }
  static MBool compare(const T& i, const MInt j) { return (i.m_partId < j); }
};

template <class T>
class sort_particleAfterCellIds {
 public:
  static MBool compare(const T& i, const T& j) { return (i.m_cellId < j.m_cellId); }
};

template <class T>
class sort_respawnParticleAfterCellIds {
 public:
  static MBool compare(const T& i, const T& j) {
    if(i.cellId == j.cellId) {
      return (i.partId < j.partId);
    }
    return (i.cellId < j.cellId);
  }
};

template <MInt nDim>
inline MBool inactiveParticle(const LPTSpherical<nDim>& particle) {
  return particle.isInvalid();
}

template <MInt nDim>
inline MBool activeParticle(const LPTSpherical<nDim>& particle) {
  return (!particle.isInvalid());
}

template <MInt nDim>
inline MBool inactiveEllipsoid(const LPTEllipsoidal<nDim>& particle) {
  return particle.isInvalid();
}

template <MInt nDim>
inline MBool activeEllipsoid(const LPTEllipsoidal<nDim>& particle) {
  return (!particle.isInvalid());
}


inline MFloat scalarProduct(const MFloat* a, const MFloat* b, const MInt length) {
  MFloat returnValue = F0;
  for(MInt count = 0; count < length; ++count) {
    returnValue += a[count] * b[count];
  }
  return returnValue;
}

inline void slerp(const MFloat* before, const MFloat* now, const MFloat time, MFloat* result, const MInt length) {
  //  TRACE();
  MFloat dotP = scalarProduct(before, now, length);
  if(dotP < 0) {
    dotP = -dotP;
    for(MInt count = 0; count < length; ++count) {
      result[count] = -now[count];
    }
  } else {
    for(MInt count = 0; count < length; ++count) {
      result[count] = now[count];
    }
  }

  if(dotP < 0.95) {
    MFloat angle = acos(dotP);
    MFloat sint1 = sin(angle * (1 - time));
    MFloat sint = sin(angle * time);
    MFloat sinA = sin(angle);
    for(MInt count = 0; count < length; ++count) {
      result[count] = (before[count] * sint1 + result[count] * sint) / sinA;
    }
  } else {
    for(MInt count = 0; count < length; ++count) {
      result[count] = before[count] * (F1 - time) + result[count] * time;
    }
  }
  maia::math::normalize(result, length);
}


/**  LPT::matrixMultiplyLeft(MFloat left[3][3], MFloat right[3][3])
 *   \brief  Matrix multiplication; matrix right is changed and contains the result
 *
 *   Feb-2011
 *   @author Rudie Kunnen
 */
inline void matrixMultiplyLeft(MFloat left[3][3], MFloat right[3][3]) {
  TRACE();

  MFloat temp[3][3];
  for(MInt i = 0; i < 3; i++) {
    for(MInt j = 0; j < 3; j++) {
      temp[i][j] = right[i][j];
      right[i][j] = 0.0;
    }
  }

  for(MInt i = 0; i < 3; i++) {
    for(MInt j = 0; j < 3; j++) {
      for(MInt k = 0; k < 3; k++) {
        right[i][j] += left[i][k] * temp[k][j];
      }
    }
  }
}

/**  LPT::matrixMultiplyRight(MFloat left[3][3], MFloat right[3][3])
 *   \brief  Matrix multiplication; matrix left is changed and contains the result
 *
 *   Feb-2011
 *   @author Rudie Kunnen
 */
inline void matrixMultiplyRight(MFloat left[3][3], MFloat right[3][3]) {
  TRACE();

  MFloat temp[3][3];
  for(MInt i = 0; i < 3; i++) {
    for(MInt j = 0; j < 3; j++) {
      temp[i][j] = left[i][j];
      left[i][j] = 0.0;
    }
  }

  for(MInt i = 0; i < 3; i++) {
    for(MInt j = 0; j < 3; j++) {
      for(MInt k = 0; k < 3; k++) {
        left[i][j] += temp[i][k] * right[k][j];
      }
    }
  }
}

/// \brief Generate a random vector in a cone defined by its opening angle
///
/// \author Sven Berger
/// \date   July 2015
/// \param[out] vec Vector of a movement in the defined cone using the provided distribution.
/// \param[in] coneAxis Center axis of the cone.
/// \param[in] length Length of the cone.
/// \param[in] openingAngle Opening angle of the cone. (Note: This means the full opening angle! It
//                          will be halfed
/// inside this function to obtain the cone angle!)
/// \param[in] PRNG random number used unknown number of times!

inline MInt randomVectorInCone(MFloat* vec, const MFloat* coneAxis, const MFloat length, const MFloat openingAngle,
                               const MInt dist, std::mt19937_64& PRNG, const MFloat distCoeff = 0.0,
                               const MFloat nozzleAngle = 0.0) {
  using namespace std;

  MFloat centerAxis[3]{};
  MFloat v1[3]{};
  MFloat v2[3]{};
  MFloat temp[3]{};
  MInt nPRNGCall = 0;

  const MFloat correctedOpeningAngle = openingAngle - 2.0 * nozzleAngle;

  std::copy_n(coneAxis, 3, &centerAxis[0]);
  // angle between center cone axis and cone shell
  MFloat distedAngle = 0;
  if(dist == PART_EMITT_DIST_GAUSSIAN) {
    //~68% of all particles are within in the center cone for distCoeff = 1.0
    normal_distribution<MFloat> distAngle(0, distCoeff * correctedOpeningAngle);
    // reject angles that are too large
    do {
      distedAngle = distAngle(PRNG);
      nPRNGCall++;
    } while(distedAngle > correctedOpeningAngle);
  } else {
    // uniform
    distedAngle = correctedOpeningAngle;
  }

  // NOTE: actually half-cone angle
  const MFloat coneAngleRad = (distedAngle / 360) * M_PI;

  ASSERT(abs(coneAxis[0]) > std::numeric_limits<MFloat>::epsilon()
             || abs(coneAxis[1]) > std::numeric_limits<MFloat>::epsilon()
             || abs(coneAxis[2]) > std::numeric_limits<MFloat>::epsilon(),
         "ERROR: ConeAxis cannot be Zero!");
  ASSERT(openingAngle <= 360, "ERROR: Opening angle cannot be larger than 360 degrees!");

  // 1. normalize centerAxis
  maia::math::normalize(&centerAxis[0], 3);

  // 2. find orthonormal basis to centerAxis
  // 2.1 choose a linear independent vector to centerAxis
  // that lies in the plane a * centerAxis = 1
  if(abs(centerAxis[0]) > 0.0) {
    v1[0] = 1.0 / centerAxis[0];
  } else if(abs(centerAxis[1]) > 0.0) {
    v1[1] = 1.0 / centerAxis[1];
  } else {
    v1[2] = 1.0 / centerAxis[2];
  }

  if(abs(coneAxis[0]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[0] = 1.0;
  } else if(abs(coneAxis[1]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[1] = 1.0;
  } else if(abs(coneAxis[2]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[2] = 1.0;
  }

  // 2.2 Use this vector to determine the first basis vector
  // by calculating v1 = v1 x centerAxis
  std::copy_n(&v1[0], 3, &temp[0]);
  maia::math::cross(&temp[0], &centerAxis[0], &v1[0]);

  // 2.3 Calculate crossproduct of v1 and centerAxis to determine the last vector v2
  maia::math::cross(&v1[0], &centerAxis[0], &v2[0]);

  // 2.4 Normalize v1 and v2
  maia::math::normalize(&v1[0], 3);
  maia::math::normalize(&v2[0], 3);

  // 3. Use orthonormal basis to determine points on spherical cap
  uniform_real_distribution<MFloat> dist0_2Pi(0, 2 * M_PI);
  // the formulation of the random number in the cos-range
  // and then taking the acos ensures the angle is in the range of:
  // -coneAngle -> coneAngle
  // when overlayed with second angle phi in the range of [ 0 -> 2 PI]
  uniform_real_distribution<MFloat> randAngle(cos(coneAngleRad), 1);

  MFloat phi = 0;
  MFloat omega = coneAngleRad;
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    // determine random point on unit circle
    phi = dist0_2Pi(PRNG);
    nPRNGCall++;

    // random angle in the given opening angle i.e in the range of [coneAngel -> 0]
    omega = acos(randAngle(PRNG));
    nPRNGCall++;
  }

  if(dist != PART_EMITT_DIST_NONE) {
    // avoid tan-formulation as omega is a random angle
    vec[0] = sin(omega) * (cos(phi) * v1[0] + sin(phi) * v2[0]) + cos(omega) * centerAxis[0];
    vec[1] = sin(omega) * (cos(phi) * v1[1] + sin(phi) * v2[1]) + cos(omega) * centerAxis[1];
    vec[2] = sin(omega) * (cos(phi) * v1[2] + sin(phi) * v2[2]) + cos(omega) * centerAxis[2];

  } else {
    // coneAngleRad must be in definition region of tan!
    // i.e. -pi/2 -> pi/2
    vec[0] = tan(omega) * (cos(phi) * v1[0] + sin(phi) * v2[0]) + centerAxis[0];
    vec[1] = tan(omega) * (cos(phi) * v1[1] + sin(phi) * v2[1]) + centerAxis[1];
    vec[2] = tan(omega) * (cos(phi) * v1[2] + sin(phi) * v2[2]) + centerAxis[2];
    // NOTE: at this point the resulting vector has the given magnitude plus the
    //      overlaying random tangential velocity component!
    //      used in some formulations of the secondary break-up!
  }

  // normalize result
  // NOTE: in this case the resulting vector has the given magnitude
  maia::math::normalize(vec, 3);


  vec[0] *= length;
  vec[1] *= length;
  vec[2] *= length;

  return nPRNGCall;
}

/// Obtain a point within a circle with a diameter around the origin.
/// \author Sven Berger
/// \date   March 2017
/// \param[out] vec Point within the defined sphere.
/// \param[in] normalDirection The circles normal direction vector.
/// \param[in] diameter Diameter of the circle.
/// \param[in] PRNG two randon numbers are generated!
inline void randomPointInCircle(MFloat* vec, const MFloat* normalDirection, const MFloat diameter,
                                std::mt19937_64& PRNG) {
  MFloat centerAxis[3]{};
  std::copy_n(normalDirection, 3, &centerAxis[0]);

  // 1. normalize centerAxis
  maia::math::normalize(&centerAxis[0], 3);

  // 2. find orthonormal basis to centerAxis
  // 2.1 choose a linear independent vector to centerAxis that lies in the plane a v1 * centerAxis = 1
  MFloat v1[3]{};
  MFloat temp[3]{};


  if(std::abs(centerAxis[0]) > 0.0) {
    v1[0] = 1.0 / centerAxis[0];
  } else if(std::abs(centerAxis[1]) > 0.0) {
    v1[1] = 1.0 / centerAxis[1];
  } else {
    v1[2] = 1.0 / centerAxis[2];
  }

  if(std::abs(normalDirection[0]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[0] = 1.0;
  } else if(std::abs(normalDirection[1]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[1] = 1.0;
  } else if(std::abs(normalDirection[2]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[2] = 1.0;
  }

  // 2.2 Use this vector to determine the first basis vector by calculating v1 = v1 x centerAxis
  std::copy_n(&v1[0], 3, &temp[0]);
  maia::math::cross(&temp[0], &centerAxis[0], &v1[0]);

  // 2.3 Calculate crossproduct of v1 and centerAxis to determine the last vector v2
  MFloat v2[3]{};
  maia::math::cross(&v1[0], &centerAxis[0], &v2[0]);

  // 2.4 Normalize v1 and v2
  maia::math::normalize(&v1[0], 3);
  maia::math::normalize(&v2[0], 3);

  // 3. Use orthonormal basis to determine points on circle

  // angle of the point within the sphere
  std::uniform_real_distribution<MFloat> dist0_2Pi(0, 2 * M_PI);
  // distance of the point from the origin within the sphere
  std::uniform_real_distribution<MFloat> randRadius(0, 1);

  MFloat phi = 0;
  MFloat radius = 0;
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    // determine random angle in unit circle
    phi = dist0_2Pi(PRNG);

    // random Radius in the unit circle
    radius = randRadius(PRNG);
  }

  vec[0] = 0.5 * diameter * radius * (cos(phi) * v1[0] + sin(phi) * v2[0]);
  vec[1] = 0.5 * diameter * radius * (cos(phi) * v1[1] + sin(phi) * v2[1]);
  vec[2] = 0.5 * diameter * radius * (cos(phi) * v1[2] + sin(phi) * v2[2]);
}

/// Obtain a random point on a circle
/// \author Sven Berger
/// \date   August 2018
/// \param[out] vec Point within the defined sphere.
/// \param[in] normalDirection The circles normal direction vector.
/// \param[in] diameter Diameter of the circle.
/// \param[in] PRNG random number used once
/// \param[in] circleSplit Split circle into this number of sections
/// \param[in] splitNo  Determine point within this partial circle
inline void randomPointOnCircle(MFloat* vec, const MFloat* normalDirection, const MFloat diameter,
                                std::mt19937_64& PRNG, const MInt circleSplit = 1, const MInt splitNo = 0) {
  MFloat centerAxis[3]{};
  std::copy_n(normalDirection, 3, &centerAxis[0]);

  // 1. normalize centerAxis
  maia::math::normalize(&centerAxis[0], 3);

  // 2. find orthonormal basis to centerAxis
  // 2.1 choose a linear independent vector to centerAxis that lies in the plane a v1 * centerAxis = 1
  MFloat v1[3]{};
  MFloat temp[3]{};


  if(std::abs(centerAxis[0]) > 0.0) {
    v1[0] = 1.0 / centerAxis[0];
  } else if(std::abs(centerAxis[1]) > 0.0) {
    v1[1] = 1.0 / centerAxis[1];
  } else {
    v1[2] = 1.0 / centerAxis[2];
  }

  if(std::abs(normalDirection[0]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[0] = 1.0;
  } else if(std::abs(normalDirection[1]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[1] = 1.0;
  } else if(std::abs(normalDirection[2]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[2] = 1.0;
  }

  // 2.2 Use this vector to determine the first basis vector by calculating v1 = v1 x centerAxis
  std::copy_n(&v1[0], 3, &temp[0]);
  maia::math::cross(&temp[0], &centerAxis[0], &v1[0]);

  // 2.3 Calculate crossproduct of v1 and centerAxis to determine the last vector v2
  MFloat v2[3]{};
  maia::math::cross(&v1[0], &centerAxis[0], &v2[0]);

  // 2.4 Normalize v1 and v2
  maia::math::normalize(&v1[0], 3);
  maia::math::normalize(&v2[0], 3);

  // 3. Use orthonormal basis to determine points on circle
  const MFloat splitCircle = (2.0 * M_PI / static_cast<MFloat>(circleSplit));

  // angle of the point within the sphere
  std::uniform_real_distribution<MFloat> dist0_2Pi(splitCircle * splitNo, splitCircle * (splitNo + 1));

  MFloat phi = 0;
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    phi = dist0_2Pi(PRNG); // determine random angle in unit circle
  }


  vec[0] = 0.5 * diameter * (cos(phi) * v1[0] + sin(phi) * v2[0]);
  vec[1] = 0.5 * diameter * (cos(phi) * v1[1] + sin(phi) * v2[1]);
  vec[2] = 0.5 * diameter * (cos(phi) * v1[2] + sin(phi) * v2[2]);
}

/// Obtain a point within a circular plane given a diameter and angle.
/// \author Sven Berger
/// \date   August 2018
/// \param[out] vec Point within the defined sphere.
/// \param[in] normalDirection The circluar plane normal direction vector.
/// \param[in] diameter Diameter of the circle.
/// \param[in] phi Angle of the point
inline void pointOnCircle(MFloat* vec, const MFloat* normalDirection, const MFloat diameter, MFloat phi) {
  MFloat centerAxis[3]{};
  std::copy_n(normalDirection, 3, &centerAxis[0]);

  // 1. normalize centerAxis
  maia::math::normalize(&centerAxis[0], 3);

  // 2. find orthonormal basis to centerAxis
  // 2.1 choose a linear independent vector to centerAxis that lies in the plane a v1 * centerAxis = 1
  MFloat v1[3]{};
  MFloat temp[3]{};


  if(std::abs(centerAxis[0]) > 0.0) {
    v1[0] = 1.0 / centerAxis[0];
  } else if(std::abs(centerAxis[1]) > 0.0) {
    v1[1] = 1.0 / centerAxis[1];
  } else {
    v1[2] = 1.0 / centerAxis[2];
  }

  if(std::abs(normalDirection[0]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[0] = 1.0;
  } else if(std::abs(normalDirection[1]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[1] = 1.0;
  } else if(std::abs(normalDirection[2]) < std::numeric_limits<MFloat>::epsilon()) {
    v1[2] = 1.0;
  }

  // 2.2 Use this vector to determine the first basis vector by calculating v1 = v1 x centerAxis
  std::copy_n(&v1[0], 3, &temp[0]);
  maia::math::cross(&temp[0], &centerAxis[0], &v1[0]);

  // 2.3 Calculate crossproduct of v1 and centerAxis to determine the last vector v2
  MFloat v2[3]{};
  maia::math::cross(&v1[0], &centerAxis[0], &v2[0]);

  // 2.4 Normalize v1 and v2
  maia::math::normalize(&v1[0], 3);
  maia::math::normalize(&v2[0], 3);

  // 3. Use orthonormal basis to determine points on circle
  vec[0] = 0.5 * diameter * (cos(phi) * v1[0] + sin(phi) * v2[0]);
  vec[1] = 0.5 * diameter * (cos(phi) * v1[1] + sin(phi) * v2[1]);
  vec[2] = 0.5 * diameter * (cos(phi) * v1[2] + sin(phi) * v2[2]);
}


///  rosin-rammler distribition function, used as initial droplet size distribution as in:
/// LARGE EDDY SIMULATION OF HIGH-VELOCITY FUEL SPRAYS:
/// STUDYING MESH RESOLUTION AND BREAKUP MODEL EFFECTS FOR SPRAY A
/// A. Wehrfritz, V. Vuorinen, O. Kaario, & M. Larmi
/// Atomization and Sprays, 23 (5): 419–442 (2013)
/// \author Tim Wegmann
/// \param[in] PRNG random number used once
/// \date   March 2021
inline MFloat rosinRammler(const MFloat min, const MFloat mean, const MFloat max, const MFloat spread,
                           std::mt19937_64& PRNG) {
  const MFloat K = 1.0 - exp(-pow((max - min) / mean, spread));

  std::uniform_real_distribution<MFloat> uni(0.0, 1.0);
  const MFloat x = uni(PRNG);

  return min + mean * pow(-log(1.0 - x * K), 1.0 / spread);
}

/// Nukiyama-Tanasawa distribution function, used as a normal-velocity ditribution for spray-wall interaction as in:
/// Spray/wall interaction models for multidimensional engine simulation
/// Z.Han, Z. Xu, N. Trigui
/// Int. J. Engine Research Vol. 1 No. 1 2000
/// \author Tim Wegmann
/// \param[in] PRNG random number used once
/// \date January 2023
inline MFloat NTDistribution(const MFloat x_mean, std::mt19937_64& PRNG) {
  std::uniform_real_distribution<MFloat> uni(0.0, 1.0);
  const MFloat x = uni(PRNG);

  return 4 / sqrt(M_PI) * POW2(x) / x_mean * exp(-POW2(x / x_mean));
}

} // namespace lpt
} // namespace maia

#endif // MAIA_PARTICLELIB_H
