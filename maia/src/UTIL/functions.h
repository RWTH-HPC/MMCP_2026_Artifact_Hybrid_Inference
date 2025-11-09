// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_FUNCTIONS_H_
#define MAIA_FUNCTIONS_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string_view>
#include <sys/stat.h>
#include <vector>
#ifdef PVPLUGIN
#include <mpi.h>
#include "src/INCLUDE/maiamacro.h"
#include "src/INCLUDE/maiatypes.h"
#include "src/compiler_config.h"
#include "src/config.h"
#else
#include "COMM/mpioverride.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "compiler_config.h"
#include "config.h"
#endif

/** \brief controlled termination of maia with a user defined error message
 *
 * \author Pascal Meysonnat, Christoph Siewert, Michael Schlottke
 * \date September/November 2011, June 2012
 */
ATTRIBUTES1(ATTRIBUTE_NORETURN)
void mTerm(const MInt errorCode, const MString& location, const MString& message = "");
#define TERM(exitval) TERMM(exitval, "")
#define TERMM(exitval, msg)                                                                                            \
  do {                                                                                                                 \
    mTerm(exitval, AT_, msg);                                                                                          \
  } while(false)

/// \brief Terminate if the given condition is fulfilled
#define TERMM_IF_COND(termCondition, msg)                                                                              \
  do {                                                                                                                 \
    if(termCondition) {                                                                                                \
      TERMM(1, msg);                                                                                                   \
    }                                                                                                                  \
  } while(false)

/// \brief Terminate if the given condition is not fulfilled
#define TERMM_IF_NOT_COND(termNotCondition, msg)                                                                       \
  do {                                                                                                                 \
    TERMM_IF_COND(!(termNotCondition), msg);                                                                           \
  } while(false)


#ifdef MAIA_TIMER_CHECKS
/// \brief Macro to be used in IO functions to assert that none of the DLB timers is running/enabled
#define CHECK_TIMERS_IO(message)                                                                                       \
  do {                                                                                                                 \
    maia::dlb::g_dlbTimerController.checkIOTimerStatus(AT_);                                                           \
  } while(false)
#else
#define CHECK_TIMERS_IO(message)                                                                                       \
  do {                                                                                                                 \
  } while(false)
#endif


/// Return the process cpu time (user time) (high-resolution timer - do not use clock())
inline MFloat cpuTime() {
  timespec tp;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp);
  const MFloat t_user = (tp.tv_sec + (MFloat)tp.tv_nsec / 1e9);
  return t_user;
}

inline MFloat wallTime() { return MPI_Wtime(); }


// Should only be used with MFloat!
template <class Real>
inline Real ABS(const Real x) {
  return std::fabs(x);
}

template <class T>
constexpr inline T mMin(const T& x, const T& y) {
  return std::min(x, y);
}
template <class T>
constexpr inline T mMax(const T& x, const T& y) {
  return std::max(x, y);
}
template <class T>
inline T MIN3(const T& x, const T& y, const T& z) {
  return std::min(std::min(x, y), z);
}
template <class T>
inline T MAX3(const T& x, const T& y, const T& z) {
  return std::max(std::max(x, y), z);
}

template <class T>
constexpr const T& Clamp(const T& v, const T& lo, const T& hi) {
  TERMM_IF_COND((hi < lo), "");
  return (v < lo) ? lo : (hi < v) ? hi : v;
}

// \brief Copies file \p fromName to file \p toName .
MInt copyFile(const MString& fromName, const MString& toName);
/// \brief Returns true if the file \p fileName exists, false otherwise.
MBool fileExists(const MString& fileName);

// ()^2
template <class Real>
constexpr Real POW2(const Real x) {
  return (x * x);
}
template <class Real>
constexpr Real POW3(const Real x) {
  return (x * x * x);
}
template <class Real>
constexpr Real POW4(const Real x) {
  return (x * x * x * x);
}
template <class Real>
constexpr Real POW5(const Real x) {
  return (x * x * x * x * x);
}
template <class Real>
constexpr Real POW6(const Real x) {
  return (x * x * x * x * x * x);
}

/** \brief      Compile time power calculation
 *  \author     Miro Gondrum
 *  \date       27.05.2021
 *  \param[in]  base
 *  \param[in]  exponent
 */
template <typename T>
T constexpr POWX(T base, MUint exponent) {
  return (exponent == 0) ? 1 : base * POWX(base, exponent - 1);
}

template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, int>::type inline constexpr signum(T x) {
  return T(0) < x;
}

template <typename T>
typename std::enable_if<std::is_signed<T>::value, int>::type inline constexpr signum(T x) {
  return (T(0) < x) - (x < T(0));
}

// the same as above but 0 -> 1
template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, int>::type inline constexpr signum0(T x) {
  return T(0) <= x;
}

template <typename T>
typename std::enable_if<std::is_signed<T>::value, int>::type inline constexpr signum0(T x) {
  return (T(0) <= x) - (x < T(0));
}


/* \brief Sutherland Law Makro for calculating the viscosity by
 *        mue = T^3/2 * (1+S/T_0)(T + S/T_0)
 *        with the default values S=110.4 k and T_0 = 273.15 K
 *
 * @author Stephan Schlimpert, 12.12.2010
 *
 * TODO labels:totest,toenhance this is WRONG and should be replaced (if T is an expression it will be
 * evaluated 3 times!)
 */
#define SUTHERLANDLAW(T) (((T)*sqrt((T)) * (m_sutherlandPlusOne)) / ((T) + (m_sutherlandConstant)))

/*  Index Macro for m_bandNghbrIdsG:
 *        IDX_LSSETDIR(cell,dir,set) =  (cell * m_maxNoSets + set ) * m_noDirs + dir
 *
 *        NOTE: ls-Solver useage only!
 *
 * @author Claudia Guenther, Jan. 2012
 */
#define IDX_LSSETDIR(i, j, k) (((i)*m_maxNoSets + (k)) * m_noDirs + (j))

/*  Index Macro for m_candidateNodeValues:
 *        IDX_LSSETNODES(cell,node,set) =  (cell * m_maxNoSets + set ) * m_noCellNodes + node
 *
 *        NOTE: fv-mb-Solver useage only!
 *
 * @author Claudia Guenther, July 2013
 */
#define IDX_LSSETNODES(i, j, k) (((i)*m_noLevelSetsUsedForMb + (k)) * m_noCellNodes + (j))


/*  Index set Macro
 *        IDX_LSSET(cell,set) =  (cell * m_maxNoSets + set )
 *
 *        NOTE: ls-Solver and cartesiansolver useage!
 *
 * @author Claudia Guenther, Jan. 2012
 */
#define IDX_LSSET(i, j) ((i)*m_maxNoSets + (j))

/*  Index Macro for m_levelSetValues:
 *        IDX_LSSETMB(cell,set) =  (cell * m_noLevelSetsUsedForMb + set )
 *   Index Macro for m_pointIsInside[ cellId ][.]:
 *        IDX_LSSETMB(point,set) =  (point * m_noLevelSetsUsedForMb + set )
 *
 *        NOTE: fv-mb-Solver useage only!
 *
 * @author Claudia Guenther, July 2013
 */
#define IDX_LSSETMB(i, j) ((i)*m_noLevelSetsUsedForMb + (j))

/// \brief The ASSERT macro is used to catch programming errors.
/// It asserts that a given condition is true at run-time
/// when compiling in debug mode.
///
/// You use it like this: ASSERT(condition,errorMsg);
/// If condition is false, the program will terminate and display errorMsg.
///
/// This is a macro for catching programming errors: i.e. errors that you
/// can make while writing code as well as errors that others can make when
/// calling your code. ASSERT can prevent them for using your code in the wrong
/// way. It can also help them understand what was wrong and how to fix it.
///
/// Do NOT use ASSERT for catching user errors at run-time!
///
/// Developer Notes:
///  - the do { } while (false) trick forces you to write ASSERT();
///  with a semicolon, this prevents evil macro things from happening.
///  - if NDEBUG is not defined if(0) will never be executed (every compiler
///  we should be using will delete the code __AFTER__ checking that the code
///  is right.
///
#ifdef MAIA_ASSERTS
#define ASSERT(condition, message)                                                                                     \
  do {                                                                                                                 \
    if(!(condition)) {                                                                                                 \
      std::cerr << "Assertion `" #condition "` failed in " << __FILE__ << " line " << __LINE__ << ": " << message      \
                << std::endl;                                                                                          \
      mTerm(1, AT_, "ASSERTION FAILED");                                                                               \
    }                                                                                                                  \
  } while(false)
#else
#define ASSERT(condition, message)                                                                                     \
  do {                                                                                                                 \
  } while(false && (condition))
#endif

#ifdef MAIA_ASSERT_ACCESSORS
#define ASSERT_ACCESSOR(condition, message) ASSERT(condition, message)
#else
#define ASSERT_ACCESSOR(condition, message)                                                                            \
  do {                                                                                                                 \
  } while(false && (condition))
#endif

namespace detail_ {
template <class T, class U>
struct APPROX_ERROR {};
} // namespace detail_

template <class T, class U>
MBool approx(const T&, const U&, const T) {
  typedef typename detail_::APPROX_ERROR<T, U>::ERROR_BOTH_TYPES_MUST_BE_MFloats error;
  error();
  return true;
}

template <>
inline MBool approx<MFloat, MFloat>(const MFloat& a, const MFloat& b, const MFloat eps) {
  return std::fabs(a - b) < eps;
}

template <class T>
MBool isApproxInt(const T&, const T) {
  typedef typename detail_::APPROX_ERROR<T, T>::ERROR_VALUE_MUST_BE_MFloat error;
  error();
  return true;
}

/// \brief Return true if argument is approximately an integer.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-08
///
/// \param[in] a Argument to check.
/// \param[in] eps Epsilon for floating point comparison.
///
/// \return True if argument is approximately an integer.
template <>
inline MBool isApproxInt<MFloat>(const MFloat& a, const MFloat eps) {
  return approx(std::fabs(std::round(a) - a), 0.0, eps);
}

/**
 * \brief Integer exponent function for non-negative exponents.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] base The base (positive or negative).
 * \param[in] exp The exponent (non-negative).
 *
 * \return base^exp, i.e. base to the power of exp
 *
 * \details Taken from http://stackoverflow.com/a/101613/1329844.
 */
inline MInt ipow(MInt base, MInt exp) {
  ASSERT(exp >= 0, "Error in ipow: exponent is less than zero (must be greater or equal to zero)!");

  MInt result = 1;
  while(exp) {
    if(exp & 1) {
#ifndef NDEBUG
      // This piece of code checks (if no -DNDEBUG flag was given during compilation) if signed integer overflow
      // occurs, and if yes, aborts the program. Signed integer overflow is detected if the result of the last
      // multiplication can not be reversed.
      MInt result_test = result * base;
      ASSERT(((base != 0 && result_test / base == result) || base == 0), "Error in ipow: signed integer overflow!");
#endif
      result *= base;
    }
    exp >>= 1;
    base *= base;
  }

  return result;
}


void checkMultiSolverGridExtents(const MInt nDim, const MFloat* centerOfGravity, const MFloat lengthLevel0,
                                 const MInt minLevel, const MFloat* targetGridCenterOfGravity,
                                 const MFloat targetGridLengthLevel0, const MInt targetGridMinLevel);

MInt loadPointCoordinatesFromFile(const MString inputFileName, const MInt nDim, std::vector<MFloat>& coordinates);

#ifndef PVPLUGIN
void writeMemoryStatistics(const MPI_Comm comm, const MInt noDomains, const MInt domainId, const MString at,
                           const MString comment = "");
#endif

inline int intSwap(int f) { // Change endian of int
  union {
    int f;
    unsigned char b[4];
  } dat1{}, dat2{};

  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}

inline float floatSwap(float f) { // Change endian of float
  union {
    float f;
    unsigned char b[4];
  } dat1{}, dat2{};

  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}

inline double doubleSwap(double f) { // Change endian of double
  union {
    double f;
    unsigned char b[8];
  } dat1{}, dat2{};

  dat1.f = f;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.f;
}

// todo labels:toenhance c++17
// inline std::string_view ltrim(std::string_view s)
//{
//    s.remove_prefix(std::distance(s.cbegin(), std::find_if(s.cbegin(), s.cend(),
//         [](int c) {return !std::isspace(c);})));
//
//    return s;
//}
//
// inline std::string_view rtrim(std::string_view s)
//{
//    s.remove_suffix(std::distance(s.crbegin(), std::find_if(s.crbegin(), s.crend(),
//        [](int c) {return !std::isspace(c);})));
//
//    return s;
//}
//
// inline std::string_view trim(std::string_view s)
//{
//    return ltrim(rtrim(s));
//}

// trim from start
static inline std::string ltrim(std::string& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) { return !std::isspace(ch); }));
  return s;
}

// trim from end
static inline std::string rtrim(std::string& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) { return !std::isspace(ch); }).base(), s.end());
  return s;
}

// trim from both ends
static inline std::string trim(std::string& s) {
  ltrim(s);
  rtrim(s);
  return s;
}

template <class ContainerT>
inline void tokenize(const std::string& str, ContainerT& tokens, const std::string& delimiters = " ",
                     MBool trimEmpty = false) {
  std::string::size_type pos = 0;
  std::string::size_type lastPos = 0;
  std::string::size_type length = str.length();

  using value_type = typename ContainerT::value_type;
  using size_type = typename ContainerT::size_type;

  while(lastPos < length + 1) {
    pos = str.find_first_of(delimiters, lastPos);
    if(pos == std::string::npos) {
      pos = length;
    }

    if(pos != lastPos || !trimEmpty) {
      tokens.push_back(value_type(str.data() + lastPos, (size_type)pos - lastPos));
    }

    lastPos = pos + 1;
  }
}

inline void createDir(const std::string& dir) {
  struct stat s {};
  if(stat(dir.c_str(), &s) < 0) {
    // Create directory if it does not yet exist
    mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXGRP | S_IXGRP);
  }
}

#endif
