// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2019 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier:    LGPL-3.0-only


#ifndef MATH_H_
#define MATH_H_

// TODO: Remove all `useEigenOld` stuff after switching to Eigen3.4.0
//#define useEigenOld

#include <algorithm>
#include <array>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <complex>
#include <functional>
#include <numeric>
#include <type_traits>
#include "COMM/mpioverride.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "debug.h"
#include "tensor.h"
namespace maia {
namespace math {
inline void quickSortImpl(MInt* a, MInt start, MInt end);
inline void quickSort(MInt* a, MInt start, MInt end);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Linear Algebra
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Adds arbritary number of vectors of dimension nDim and returns result in vector result. Vector
/// result can contain an offset. Usage:
/// Say you want to add vectors a, b, c of dimension 3 and store result in a: vecAdd<3>(a, b, c)
template <MInt nDim, typename T, typename U>
inline T* vecAdd(T* const result, const U* const a) {
  std::transform(a, a + nDim, result, result, std::plus<T>());
  return result;
}
template <MInt nDim, typename T, typename U, typename... Ts>
inline T* vecAdd(T* const result, const U* const a, const Ts... b) {
  std::transform(a, a + nDim, result, vecAdd<nDim>(result, b...), std::plus<T>());
  return result;
}

/// Computes the arithmetic mean of an arbritary number of vectors of dimension nDim and returns it in M
template <MInt nDim, typename T, typename... Ts>
inline void vecAvg(T* const M, const Ts* const... args) {
  std::fill_n(M, nDim, 0.0);
  maia::math::vecAdd<nDim>(M, args...);
  constexpr std::size_t n = sizeof...(Ts);
  std::transform(M, M + nDim, M, std::bind(std::multiplies<T>(), std::placeholders::_1, 1.0 / n));
}


/// Cross product C = U x V [3D]
/// \tparam T Type
/// \param[in] u U
/// \param[in] v V
/// \param[out] c C
/// \author Sven Berger
template <typename T>
inline void cross(const T* const u, const T* const v, T* const c) {
  c[0] = u[1] * v[2] - v[1] * u[2];
  c[1] = u[2] * v[0] - v[2] * u[0];
  c[2] = u[0] * v[1] - v[0] * u[1];
}

/// Cross product C = U x V [3D]
/// \tparam T Type
/// \param[in] u U
/// \param[in] v V
/// \return New copy of an array
template <typename T>
inline std::array<T, 3> cross(const std::array<T, 3>& u, const std::array<T, 3>& v) {
  std::array<T, 3> result;
  cross(&u[0], &v[0], &result[0]);
  return result;
}
// 2D case
template <typename T>
inline T cross(const std::array<T, 2>& u, const std::array<T, 2>& v) {
  return u[0] * v[1] - u[1] * v[0];
}

/// Cross product C = U x V [3D]
/// \tparam T Type
/// \param[in] u U
/// \param[in] v V
/// \return New copy of an array
template <typename T>
inline T* cross(const T (&u)[3], const T (&v)[3]) {
  T result[3]{};
  cross(&u[0], &v[0], &result[0]);
  return result;
}
// 2D case
template <typename T>
inline T cross(const T (&u)[2], const T (&v)[2]) {
  return u[0] * v[1] - u[1] * v[0];
}

/// Length of a vector (L2-norm)
/// \tparam T Type
/// \tparam N Length of vector
/// \param u Vector
/// \return Magnitude of vector
/// \author Sven Berger
template <typename T, std::size_t N>
inline MFloat norm(const std::array<T, N>& u) {
  static_assert(N > 1, "ERROR: Invalid norm call!");

  if(N == 2) {
    return std::sqrt(u[0] * u[0] + u[1] * u[1]);
  }
  if(N == 3) {
    return std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  }
  return std::sqrt(std::accumulate(u.begin(), u.end(), 0.0, [](const MFloat& a, const T& b) { return a + b * b; }));
}

/// Length of a vector (L2-norm)
/// \tparam T Type
/// \param N Length of vector
/// \param u Vector
/// \return Magnitude of vector
/// \author Sven Berger
template <typename T>
inline MFloat norm(const T* const u, const MInt N) {
  ASSERT(N > 1, "ERROR: Invalid norm call!");

  if(N == 2) {
    return std::sqrt(u[0] * u[0] + u[1] * u[1]);
  }
  if(N == 3) {
    return std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  }

  MFloat tmp = 0;
  for(MInt i = 0; i < N; i++) {
    tmp += u[i] * u[i];
  }

  return std::sqrt(tmp);
}

///  Normalize a given vector (inplace)
/// \tparam T Type
/// \tparam N Length of the vector
/// \param u Vector to be normalized
/// \author Sven Berger
template <typename T, std::size_t N>
inline void normalize(std::array<T, N>& u) {
  const MFloat inverse = 1.0 / norm(u);

  if(N == 2) {
    u[0] *= inverse;
    u[1] *= inverse;
    return;
  }
  if(N == 3) {
    u[0] *= inverse;
    u[1] *= inverse;
    u[2] *= inverse;
    return;
  }

  std::for_each(u.begin(), u.end(), [inverse](T& a) { a *= inverse; });
}

///  Normalize a given vector (inplace)
/// \tparam T Type
/// \param N Length of the vector
/// \param u Vector to be normalized
/// \author Sven Berger
template <typename T>
inline void normalize(T* const u, const MInt N) {
  const MFloat inverse = 1.0 / norm(u, N);

  if(N == 2) {
    u[0] *= inverse;
    u[1] *= inverse;
    return;
  }
  if(N == 3) {
    u[0] *= inverse;
    u[1] *= inverse;
    u[2] *= inverse;
    return;
  }

  for(MInt i = 0; i < N; i++) {
    u[i] *= inverse;
  }
}

///  Normalize a given vector
/// \tparam T Type
/// \tparam N Length of the vector
/// \param u Vector to be normalized
/// \return New normalized vector
/// \author Sven Berger
template <typename T, std::size_t N>
inline std::array<T, N> normalized(const std::array<T, N>& u) {
  std::array<T, N> n(u);
  normalize(n);
  return n;
}

template <MInt nDim>
inline MFloat distance(const MFloat* a, const MFloat* b) {
  IF_CONSTEXPR(nDim == 2) { return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1])); }
  IF_CONSTEXPR(nDim == 3) { return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1]) + POW2(a[2] - b[2])); }
  return -NAN;
}

inline MFloat distance(const std::array<MFloat, 2> a, const std::array<MFloat, 2> b) {
  return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1]));
}

inline MFloat distance(const std::array<MFloat, 3> a, const std::array<MFloat, 3> b) {
  return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1]) + POW2(a[2] - b[2]));
}

inline MFloat distance(const MFloat* const a, const std::array<MFloat, 2> b) {
  return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1]));
}

inline MFloat distance(const MFloat* const a, const std::array<MFloat, 3> b) {
  return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1]) + POW2(a[2] - b[2]));
}

inline MFloat distance(const std::array<MFloat, 2> a, const MFloat* const b) {
  return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1]));
}

inline MFloat distance(std::array<MFloat, 3> a, const MFloat* const b) {
  return sqrt(POW2(a[0] - b[0]) + POW2(a[1] - b[1]) + POW2(a[2] - b[2]));
}

inline void multiplyMatricesSq(MFloatScratchSpace& m1, MFloatScratchSpace& m2, MFloatScratchSpace& result, MInt dim) {
  // init
  for(MInt i = 0; i < dim; i++) {
    for(MInt j = 0; j < dim; j++) {
      result(i, j) = 0.0;
    }
  }

  // multiply
  for(MInt i = 0; i < dim; i++) {
    for(MInt j = 0; j < dim; j++) {
      for(MInt k = 0; k < dim; k++) {
        result(i, j) += m1(i, k) * m2(k, j);
      }
    }
  }
}

inline void multiplyMatrices(MFloatScratchSpace& m1, MFloatScratchSpace& m2, MFloatScratchSpace& result, MInt m1_n,
                             MInt m1_m, MInt m2_n, MInt m2_m) {
  // m1 has the dimension n x m
  // m2 has the dimension n x m
  if(m1_m != m2_n) mTerm(1, AT_, "Dimension of matrix multiplication does not match!!!");
  // init
  for(MInt i = 0; i < m1_n; i++) {
    for(MInt j = 0; j < m2_m; j++) {
      result(i, j) = F0;
    }
  }

  // multiply
  for(MInt i = 0; i < m1_n; i++) {
    for(MInt j = 0; j < m2_m; j++) {
      for(MInt k = 0; k < m1_m; k++) {
        result(i, j) += m1(i, k) * m2(k, j);
      }
    }
  }
}

inline void addMatrices(MFloatScratchSpace& m1, MFloatScratchSpace& m2, MFloatScratchSpace& result, MInt dim1,
                        MInt dim2) {
  for(MInt i = 0; i < dim1; i++) {
    for(MInt j = 0; j < dim2; j++) {
      result(i, j) = m1(i, j) + m2(i, j);
    }
  }
}

inline MFloat frobeniusMatrixNormSquared(MFloatScratchSpace& m, MInt dim1, MInt dim2) {
  MFloat ret_val = 0.0;
  for(MInt i = 0; i < dim1; i++) {
    for(MInt j = 0; j < dim2; j++) {
      ret_val += m(i, j) * m(i, j);
    }
  }

  return ret_val;
}

inline MFloat frobeniusMatrixNorm(MFloatScratchSpace& m, MInt dim1, MInt dim2) {
  return sqrt(frobeniusMatrixNormSquared(m, dim1, dim2));
}

inline MFloat besselJ0(MFloat x) {
  // This subroutine calculates the First Kind Bessel Function of
  // order 0, for any real number X. The polynomial approximation by
  // series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
  // REFERENCES:
  // M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
  // C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
  // VOL.5, 1962.

  const MFloat p1 = 1.0, p2 = -0.1098628627e-2, p3 = 0.2734510407e-4, p4 = -0.2073370639e-5, p5 = 0.2093887211e-6,
               q1 = -0.1562499995e-1, q2 = 0.1430488765e-3, q3 = -0.6911147651e-5, q4 = 0.7621095161e-6,
               q5 = -0.9349451520e-7, r1 = 57568490574.0, r2 = -13362590354.0, r3 = 651619640.7, r4 = -11214424.18,
               r5 = 77392.33017, r6 = -184.9052456, s1 = 57568490411.0, s2 = 1029532985.0, s3 = 9494680.718,
               s4 = 59272.64853, s5 = 267.8532712, s6 = 1.0;
  MFloat ax, fr, fs, z, fp, fq1, xx, y, tmp;

  if(approx(x, 0.0, std::numeric_limits<MFloat>::epsilon())) {
    return 1.0;
  }

  ax = fabs(x);
  if(ax < 8.0) {
    y = x * x;
    fr = r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6))));
    fs = s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6))));
    tmp = fr / fs;
  } else {
    z = 8. / ax;
    y = z * z;
    xx = ax - 0.785398164;
    fp = p1 + y * (p2 + y * (p3 + y * (p4 + y * p5)));
    fq1 = q1 + y * (q2 + y * (q3 + y * (q4 + y * q5)));
    tmp = std::sqrt(0.636619772 / ax) * (fp * cos(xx) - z * fq1 * sin(xx));
  }
  return tmp;
}

inline MFloat besselJ1(MFloat x) {
  // This subroutine calculates the First Kind Bessel Function of
  // order 1, for any real number X. The polynomial approximation by
  // series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
  // REFERENCES:
  // M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
  // C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
  // VOL.5, 1962.

  const MFloat p1 = 1.0, p2 = 0.183105e-2, p3 = -0.3516396496e-4, p4 = 0.2457520174e-5, p5 = -0.240337019e-6,
               p6 = 0.636619772, q1 = 0.04687499995, q2 = -0.2002690873e-3, q3 = 0.8449199096e-5, q4 = -0.88228987e-6,
               q5 = 0.105787412e-6, r1 = 72362614232.0, r2 = -7895059235.0, r3 = 242396853.1, r4 = -2972611.439,
               r5 = 15704.48260, r6 = -30.16036606, s1 = 144725228442.0, s2 = 2300535178.0, s3 = 18583304.74,
               s4 = 99447.43394, s5 = 376.9991397, s6 = 1.0;

  MFloat ax, fr, fs, y, z, fp, fq1, xx, tmp;

  ax = fabs(x);
  if(ax < 8.0) {
    y = x * x;
    fr = r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6))));
    fs = s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6))));
    tmp = x * (fr / fs);
  } else {
    z = 8.0 / ax;
    y = z * z;
    xx = ax - 2.35619491;
    fp = p1 + y * (p2 + y * (p3 + y * (p4 + y * p5)));
    fq1 = q1 + y * (q2 + y * (q3 + y * (q4 + y * q5)));
    tmp = sqrt(p6 / ax) * (cos(xx) * fp - z * sin(xx) * fq1) * (x < 0.0 ? -fabs(s6) : fabs(s6));
  }
  return tmp;
}

// bessel functions of first kind
inline MFloat besselJi(MInt order, MFloat x) {
  if(order == 0) {
    return besselJ0(x);
  }
  if(order == 1) {
    return besselJ1(x);
  }
  mTerm(1, AT_, "invalid order");
}

// cosinus using a lookuptable
inline MFloat lincos(MFloat arg) {
  MFloat SIGN = F1;
  // lookuptable has 91 entries from 0 to pi/2
  const MFloat darg = PIB2 / 90;

  arg = std::fmod(std::abs(arg), F2 * PI);

  if(arg > PI) {
    SIGN *= -F1;
    arg -= PI;
  }

  if(arg > PIB2) {
    SIGN *= -F1;
    arg = PI - arg;
  }

  const auto lowerInd = static_cast<MInt>(arg / darg);
  const MFloat lowerVal = FTRIG[lowerInd];
  const MFloat higherVal = FTRIG[lowerInd + 1];
  const MFloat deltaArg = arg - lowerInd * darg;
  const MFloat fac = deltaArg / darg;

  return (lowerVal * (F1 - fac) + higherVal * fac) * SIGN;
}

/** void sortEigenVectors(MFloat A[3][3], MFloat w[3])
 *
 *  \brief Sorts the eigenvalues and the associated eigenvectors from large to small
 *
 *         Parameters:
 *         A: The symmetric input matrix
 *         w: Storage buffer for eigenvalues
 */
inline void sortEigenVectors(MFloat A[3][3], MFloat w[3]) {
  const MInt dim = 3;
  MInt k;
  for(MInt i = 0; i < (dim - 1); ++i) {
    MFloat p = w[k = i];
    for(MInt j = i; j < dim; ++j) {
      if(w[j] >= p) {
        p = w[k = j];
      }
    }
    if(k != i) {
      w[k] = w[i];
      w[i] = p;
      for(MInt j = 0; j < dim; ++j) {
        p = A[j][i];
        A[j][i] = A[j][k];
        A[j][k] = p;
      }
    }
  }
}

inline void quatMult(const MFloat* const qA, const MFloat* const qB, MFloat* const qC) {
  MFloat crossProduct[3]{};
  cross(qA, qB, &crossProduct[0]);

  const MFloat dot = std::inner_product(qA, &qA[3], qB, 0.0);

  for(MInt n = 0; n < 3; n++) {
    qC[n] = qA[3] * qB[n] + qB[3] * qA[n] + crossProduct[n];
  }
  qC[3] = qA[3] * qB[3] - dot;
}

template <typename T>
inline MInt sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/** \brief help-function for engine calculations which returns the crank-angle
 *         for a given time, Strouhals-number and crank-angle offset
 *         mode = 0: return CAD in range of (0-720)
 *         mode = 1: return accumulated crankAnge in radian
 *
 *  \author Tim Wegmann
 */
inline MFloat crankAngle(const MFloat time, const MFloat Strouhal, const MFloat offset, const MInt mode) {
  const MFloat mu2 = Strouhal * F2 * PI;
  MFloat cad = mu2 * time;

  if(mode == 0) {
    const MInt maxNoCycles = 20;
    for(MInt cycle = maxNoCycles; cycle > 0; cycle--) {
      if(cad >= 4 * PI * cycle) {
        cad = cad - 4 * PI * cycle;
      }
    }
    cad = cad * 180 / PI;

    cad = cad + offset;

  } else {
    ASSERT(mode == 1, "Incorrect mode!");
    cad = cad + offset * PI / 180;
  }

  return cad;
}

/**
 * \brief rotation matrix co-rotating(~inertial) frame -> body-fixed frame
 * \author Lennart Schneiders
 */
inline void computeRotationMatrix(MFloatScratchSpace& R, MFloat* q) {
  R(0, 0) = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  R(1, 1) = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  R(2, 2) = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];

  R(0, 1) = F2 * (q[1] * q[2] + q[0] * q[3]);
  R(0, 2) = F2 * (q[1] * q[3] - q[0] * q[2]);
  R(1, 0) = F2 * (q[1] * q[2] - q[0] * q[3]);
  R(1, 2) = F2 * (q[2] * q[3] + q[0] * q[1]);
  R(2, 0) = F2 * (q[1] * q[3] + q[0] * q[2]);
  R(2, 1) = F2 * (q[2] * q[3] - q[0] * q[1]);
}

/**
 * \brief
 * \author
 */
inline void rotation2quaternion(MFloat* rotation, MFloat* quaternion) {
  quaternion[0] = cos(F1B2 * rotation[1]) * cos(F1B2 * (rotation[2] + rotation[0]));
  quaternion[1] = sin(F1B2 * rotation[1]) * sin(F1B2 * (rotation[2] - rotation[0]));
  quaternion[2] = sin(F1B2 * rotation[1]) * cos(F1B2 * (rotation[2] - rotation[0]));
  quaternion[3] = cos(F1B2 * rotation[1]) * sin(F1B2 * (rotation[2] + rotation[0]));
}

/**
 * \brief c=A*b
 * \author Lennart Schneiders
 */
inline void matrixVectorProduct(MFloat* c, MFloatScratchSpace& A, MFloat* b) {
  for(MInt i = 0; i < A.size0(); i++) {
    c[i] = F0;
    for(MInt j = 0; j < A.size1(); j++) {
      c[i] += A(i, j) * b[j];
    }
  }
}

/**
 * \brief c=A^t*b
 * \author Lennart Schneiders
 */
inline void matrixVectorProductTranspose(MFloat* c, MFloatScratchSpace& A, MFloat* b) {
  for(MInt i = 0; i < A.size1(); i++) {
    c[i] = F0;
    for(MInt j = 0; j < A.size0(); j++) {
      c[i] += A(j, i) * b[j];
    }
  }
}


// -----------------------------------------------------------------------------------------------------------


inline MInt inverse(MFloat** a, MFloat** ainv, MInt n, const MFloat epsilon) {
  MInt s;
  MInt pRow = 0; // pivot row
  MBool error = false;
  MFloat f;
  MFloat maximum;
  MInt pivot = 1;

  // add the unity matrix
  for(MInt i = 0; i < n; i++) {
    for(MInt j = 0; j < n; j++) {
      a[i][n + j] = F0;
      if(i == j) {
        a[i][n + j] = F1;
      }
    }
  }

  // Gauss algorithm
  error = false;
  s = 0;
  while(s < n) {
    maximum = fabs(a[s][s]);
    if(pivot) {
      pRow = s;
      for(MInt i = s + 1; i < n; i++) {
        if(fabs(a[i][s]) > maximum) {
          maximum = fabs(a[i][s]);
          pRow = i;
        }
      }
    }
    if(maximum < epsilon) {
      error = true;
    }

    if(error) {
      std::cerr << "Error in matrix inverse computation " << s << " " << a[s][s] << std::endl;
      for(MInt i = 0; i < n; i++) {
        for(MInt j = 0; j < n; j++) {
          std::cerr << a[i][j] << " ";
        }
        std::cerr << std::endl;
      }
      // mTerm(1, AT_, "Error in matrix inverse computation ");
      std::cerr << "Error in matrix inverse computation " << std::endl;
      return 0;
    }

    if(pivot) {
      if(pRow != s) // exchange rows if required
      {
        MFloat h;
        for(MInt j = s; j < 2 * n; j++) {
          h = a[s][j];
          a[s][j] = a[pRow][j];
          a[pRow][j] = h;
        }
      }
    }

    f = a[s][s];
    for(MInt j = s; j < 2 * n; j++) {
      a[s][j] = a[s][j] / f;
    }

    // elimination
    for(MInt i = 0; i < n; i++) {
      if(i != s) {
        f = -a[i][s];
        for(MInt j = s; j < 2 * n; j++) {
          a[i][j] += f * a[s][j];
        }
      }
    }
    s++;
  }

  if(error) {
    std::cerr << "Error 2 in inverse matrix computation" << std::endl;
    // mTerm(1,AT_,"Error 2 in inverse matrix computation");
    return 0;
  }

  // copy
  for(MInt i = 0; i < n; i++) {
    for(MInt j = 0; j < n; j++) {
      ainv[i][j] = a[i][n + j];
    }
  }

  return 1;
}

inline MInt quickSortPartition(MInt* a, MInt start, MInt end) {
  MInt i, temp;
  //---

  i = start - 1;

  for(MInt j = start; j < end; j++) {
    if(a[j] <= a[end]) {
      i++;
      // swap elements i and j
      temp = a[i];
      a[i] = a[j];
      a[j] = temp;
    }
  }
  // swap the pivot
  temp = a[i + 1];
  a[i + 1] = a[end];
  a[end] = temp;

  return i + 1;
}

inline void quickSortImpl(MInt* a, MInt start, MInt end) {
  MInt pivot;

  if(start >= end) {
    return;
  }

  pivot = maia::math::quickSortPartition(a, start, end);
  maia::math::quickSort(a, start, pivot - 1);
  maia::math::quickSort(a, pivot + 1, end);
}

inline void quickSort(MInt* a, MInt start, MInt end) { maia::math::quickSortImpl(a, start, end); }


// removes double entries in arrays
// retains the sorting
// returns the new size of the array
inline MInt removeDoubleEntries(MInt* a, MInt size) {
  for(MInt i = 1; i < size; i++) {
    // find double entry
    if(a[i] != a[i - 1]) {
      continue;
    }

    // shift entries forward
    for(MInt k = i; k < size - 1; k++) {
      a[k] = a[k + 1];
    }
    size--;
    i--;
  }

  return size;
}

/**
 * \brief Evaluates the Legendre polynomial and its derivative of degree Nmax at
 *        point x.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-06
 *
 * \param[in] Nmax The polynomial degree.
 * \param[in] x The evaluation point.
 * \param[out] polynomial The resulting value of the Legendre polynomial.
 * \param[out] derivative The resulting value of the derivative.
 *
 * \details Taken from Kopriva09, p. 63, algorithm 22.
 */
inline void calculateLegendrePolyAndDeriv(MInt Nmax, MFloat x, MFloat* polynomial, MFloat* derivative) {
  // TRACE();

  // Create temporary storage for the polynomial and its derivative
  MFloat poly = F0;
  MFloat deriv = F0;

  // Explicitly calculate the first two values of the three term recursion
  if(Nmax == 0) {
    poly = F1;
    deriv = F0;
  } else if(Nmax == 1) {
    poly = x;
    deriv = F1;
  } else {
    MFloat polyLast1, polyLast2, derivLast1, derivLast2;

    polyLast2 = F1;
    polyLast1 = x;
    derivLast2 = F0;
    derivLast1 = F1;

    // Calculate the polynomial and its derivative for higher degrees
    // (Nmax >= 2)
    for(MInt k = 2; k <= Nmax; k++) {
      poly = (F2 * k - F1) / k * x * polyLast1 - (k - F1) / k * polyLast2;
      deriv = derivLast2 + (F2 * k - F1) * polyLast1;
      polyLast2 = polyLast1;
      polyLast1 = poly;
      derivLast2 = derivLast1;
      derivLast1 = deriv;
    }
  }

  // Save results to pointer locations
  *polynomial = poly;
  *derivative = deriv;
}

/**
 * \brief Calculate the Gauss integration nodes and weight for the Legendre
 *        polynomials on the interval [-1,1].
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] Nmax Maximum degree of the polynomials.
 * \param[out] nodes The resulting integration nodes.
 * \param[out] wInt The resulting integration weights.
 *
 * \details Taken from Kopriva09, p. 64, algorithm 23.
 */
inline void calculateLegendreGaussNodesAndWeights(MInt Nmax, MFloat* nodes, MFloat* wInt) {
  TRACE();

  // Reset nodes and weights
  const MInt noNodes = Nmax + 1;
  std::fill_n(&nodes[0], noNodes, F0);
  std::fill_n(&wInt[0], noNodes, F0);

  // Set tolerance and number of iterations. According to Kopriva09,
  // 1.0E-15 and 10 should be more than sufficient.
  const MFloat tol = F4 * MFloatEps;
  const MInt noIterations = 10;

  // Catch simple cases before going into the full loop
  if(Nmax == 0) {
    nodes[0] = F0;
    wInt[0] = F2;
  } else if(Nmax == 1) {
    nodes[0] = -sqrt(F1B3);
    wInt[0] = 1;
    nodes[1] = -nodes[0];
    wInt[1] = wInt[0];
  } else {
    // Use symmetry property of the roots of the Legendre polynomials
    for(MInt j = 0; j < (Nmax + 1) / 2; j++) {
      // Calculate starting guess for Newton method
      nodes[j] = -cos((F2 * j + F1) / (F2 * Nmax + F2) * PI);

      // Use Newton method to find root of Legendre polynomial
      // -> this is also the integration node
      MFloat poly, deriv;
      for(MInt k = 0; k < noIterations; k++) {
        calculateLegendrePolyAndDeriv(Nmax + 1, nodes[j], &poly, &deriv);
        MFloat delta = -poly / deriv;
        nodes[j] += delta;

        // Stop iterations if error is small enough
        if(fabs(delta) <= tol * fabs(nodes[j])) {
          break;
        }
      }

      // Calculate weight
      calculateLegendrePolyAndDeriv(Nmax + 1, nodes[j], &poly, &deriv);
      wInt[j] = F2 / ((1 - nodes[j] * nodes[j]) * deriv * deriv);

      // Set nodes and weights according to symmetry properties
      nodes[Nmax - j] = -nodes[j];
      wInt[Nmax - j] = wInt[j];
    }
  }

  // If odd number of nodes (noNodes = Nmax + 1), set center node to
  // origin (0.0) and calculate weight
  if(Nmax % 2 == 0) {
    MFloat poly, deriv;
    calculateLegendrePolyAndDeriv(Nmax + 1, F0, &poly, &deriv);
    nodes[Nmax / 2] = F0;
    wInt[Nmax / 2] = F2 / (deriv * deriv);
  }
}
/**
 * \brief radial base function
 * \author Lennart Schneiders
 * \note important: provide rapid decrease of return value with R, otherwise stability issues with weighted least
 * squares when inlcuding diagonal neighbor cells! \note suggested value for R0 is cellLength
 */
inline MFloat RBF(const MFloat R, const MFloat R0) { return (F1 / sqrt(F1 + POW2(R / R0))); }


// ---------------------------------------------------------------------------


/**
 * \brief series of transition functions (such that deltaFun(r<=r0)=0 and deltaFun(r>=r1)=1) with ascending smoothness
 * \author Lennart Schneiders
 */
#define MAIA_TRANSITION_FUNCTION 2

inline MFloat deltaFun(const MFloat r, const MFloat r0, const MFloat r1) {
  MFloat R = mMin(F1, mMax(F0, ((r - r0) / (r1 - r0))));
#if MAIA_TRANSITION_FUNCTION == 0
  return (r < r1) ? F0 : F1;
#elif MAIA_TRANSITION_FUNCTION == 1
  return R;
#elif MAIA_TRANSITION_FUNCTION == 2
  return POW2(R) * (3.0 - 2.0 * R);
#elif MAIA_TRANSITION_FUNCTION == 3
  return POW3(R) * (10.0 + R * (6.0 * R - 15.0)); // may cause stability issues
#elif MAIA_TRANSITION_FUNCTION == 4
  return POW3(R) * (6.0 + R * (2.0 * R - 7.0)); // approximates 0.5*(1.0-cos(PI*pow(R,1.25)));
#elif MAIA_TRANSITION_FUNCTION == 5
  return F1B2 * (F1 + sin(PI * (F3B2 + R)));
#endif
}

// MATLAB linspace
inline std::vector<MFloat> linSpace(const MFloat start, const MFloat end, const MInt num) {
  std::vector<MFloat> linspaced(num);
  MFloat delta = (end - start) / (MFloat(num) - F1);
  for(auto i = 0; i < num; i++) {
    linspaced[i] = start + delta * i;
  }
  return linspaced;
}

inline MFloat getSector(MFloat y, MFloat z, MFloat azimuthalAngle) {
  MFloat angle = atan(z / y);
  if(y < 0) {
    if(z >= 0) {
      angle += PI;
    } else {
      angle -= PI;
    }
  }
  angle += PI;
  return (angle - std::fmod(angle, azimuthalAngle / 180.0 * PI)) / (azimuthalAngle / 180.0 * PI);
}

inline MFloat getAngle(MFloat y, MFloat z) {
  MFloat angle = atan(z / y);
  if(y < 0) {
    if(z >= 0) {
      angle += PI;
    } else {
      angle -= PI;
    }
  }
  return angle;
}

//------------------------------------------------------------------------

template <typename T, std::size_t N>
MFloat determinant(std::array<T, N>& m);
template <typename T, std::size_t N>
MFloat determinant(std::array<std::array<T, N>, N>& m);
void invert(MFloat* A, const MInt m, const MInt n);
template <class T>
void invert(T& A, T& AInv, const MInt m, const MInt n);
template <class T>
void invert(T& A, T& weights, T& AInv, const MInt m, const MInt n);
template <class T>
MInt invertR(T& A, T& weights, T& AInv, const MInt m, const MInt n);
void solveDenseMatrix(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x);
void solveSparseMatrixIterative(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x);
void solveSparseMatrix(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x);
void multiplySparseMatrixVector(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b_in, MFloat* x_final);
void calcEigenValues(MFloat A[3][3], MFloat w[3]);
void calcEigenValues(MFloat** A_in, MFloat* lambda_in, const MInt m);
void calcEigenVectors(MFloat A[3][3], MFloat Q[3][3], MFloat w[3]);
template <MInt nDim>
void solveQR(std::array<std::array<MFloat, nDim>, nDim>& A_, std::array<MFloat, nDim>& b_);
template <typename T, std::size_t N>
void adjointRow(std::array<std::array<T, N>, N>& m, std::array<T, N>& A, const MInt r);
template <typename T>
void adjoint1stRow4x4(std::array<std::array<T, 4>, 4>& m, std::array<T, 4>& A);
template <typename T, std::size_t N>
void adjoint1stRow(std::array<std::array<T, N>, N>& m, std::array<T, N>& A);
std::vector<MFloat> svd(MFloat* const A, MFloat* const b, const MInt m, const MInt n, MFloat* const x);
} // namespace math
} // namespace maia
#endif // MATH_H_
