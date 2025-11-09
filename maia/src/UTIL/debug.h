// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DEBUG_H
#define DEBUG_H

#include <cstdlib>
#include <iostream>
#include <sstream>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "IO/infoout.h"
#include "UTIL/timer.h"
#include "globalvariables.h"

#define __FUNCTION_LOCATION__ std::string(__FILE__) + ": " + std::string(__FUNCTION__)

#ifdef MAIA_PROFILING
// Define profiling macros
#define DEBUG_DISPLAY_ON
#define DEBUG_DISPLAY_OFF
#define SET_DEBUG_LEVEL(a) MDebug::setLevelNotEnabled(a)
#define DEBUG(a, b)
#define PROFILE(id, a, b)                                                                                              \
  do {                                                                                                                 \
    if(b == MAIA_DEBUG_TRACE_IN) {                                                                                     \
      Profile::s_functionTimings[id].in();                                                                             \
    } else if(b == MAIA_DEBUG_TRACE_OUT) {                                                                             \
      Profile::s_functionTimings[id].out();                                                                            \
    }                                                                                                                  \
  } while(false)

// Note: this needs to be without braces/do-while, using this will cause the Tracer object to directly go out of of
// scope again
// Note: to enable profiling of all TRACE-instrumented functions set the last argument to 'true'
#define TRACE()                                                                                                        \
  static MInt timingId = Profile::getTimingId(FUN_);                                                                   \
  maia::debug::Tracer USE_ONLY_ONCE_PER_METHOD(FUN_, LOC_, timingId, false)
// Use TRACE_PROFILE instead of TRACE if only certain functions should be profiled
#define TRACE_PROFILE()                                                                                                \
  static MInt timingId = Profile::getTimingId(FUN_);                                                                   \
  maia::debug::Tracer USE_ONLY_ONCE_PER_METHOD(FUN_, LOC_, timingId, true)

#elif defined MAIA_DEBUG_FUNCTION
// Define debugging macros
#define DEBUG_DISPLAY_ON MDebug::displayOn();
#define DEBUG_DISPLAY_OFF MDebug::displayOff();
#define SET_DEBUG_LEVEL(a) MDebug::setLevel(a);
#define DEBUG(a, b)                                                                                                    \
  do {                                                                                                                 \
    std::ostringstream message;                                                                                        \
    message << a;                                                                                                      \
    MDebug::maiaprint(message.str(), b);                                                                               \
  } while(false)
#define PROFILE(id, a, b)
#define TRACE() maia::debug::Tracer USE_ONLY_ONCE_PER_METHOD(FUN_, LOC_)
#define TRACE_PROFILE() TRACE()

#else
// Reset all profiling/debugging macros
#define DEBUG_DISPLAY_ON
#define DEBUG_DISPLAY_OFF
#define SET_DEBUG_LEVEL(a) MDebug::setLevelNotEnabled(a)
#define DEBUG(a, b)
#define PROFILE(id, a, b)
#define TRACE()
#define TRACE_PROFILE()
#endif


// Additional debugging macros
#ifdef MAIA_DEBUG_FUNCTION
#define DOUT(a)                                                                                                        \
  do {                                                                                                                 \
    std::cerr << "DEBUG: VALUE: '" << #a << "' = '" << a << "' at " << AT_ << std::endl;                               \
  } while(false)

#define DOUT_IF(condition, a)                                                                                          \
  do {                                                                                                                 \
    if(condition) {                                                                                                    \
      std::cerr << "DEBUG: VALUE: '" << #a << "' = '" << a << "' at " << AT_ << std::endl;                             \
    }                                                                                                                  \
  } while(false)

#define DLOC()                                                                                                         \
  do {                                                                                                                 \
    std::cerr << "DEBUG: DOMAIN: " << globalDomainId() << "; LOCATION: " << AT_ << std::endl;                          \
  } while(false)
#define DMSG(a)                                                                                                        \
  do {                                                                                                                 \
    std::cerr << "DEBUG: " << globalDomainId() << "; MESSAGE: '" << a << "' at " << AT_ << std::endl;                  \
  } while(false)
#else
// Define empty macros with do-while clause to enfore a semicolon after the
// macro
#define DOUT(a)                                                                                                        \
  do {                                                                                                                 \
  } while(false)
#define DOUT_IF(condition, a)                                                                                          \
  do {                                                                                                                 \
  } while(false)
#define DLOC()                                                                                                         \
  do {                                                                                                                 \
  } while(false)
#define DMSG(a)                                                                                                        \
  do {                                                                                                                 \
  } while(false)
#endif

//! This enum holds the error levels of the DEBUG class.
using MDebugLevel = enum {
  MAIA_DEBUG_ASSERTION = 1,
  MAIA_DEBUG_IO = 2,
  MAIA_DEBUG_ALLOCATION = 4,
  MAIA_DEBUG_TRACE = 8,
  MAIA_DEBUG_TRACE_IN = 16,
  MAIA_DEBUG_TRACE_OUT = 32,
  MAIA_DEBUG_LEVEL1 = 64,
  MAIA_DEBUG_LEVEL2 = 128,
  MAIA_DEBUG_LEVEL3 = 256,
  MAIA_DEBUG_LEVEL4 = 512,
  MAIA_DEBUG_LEVEL5 = 1024,
  MAIA_DEBUG_USER1 = 2048,
  MAIA_DEBUG_USER2 = 4096,
  MAIA_DEBUG_USER3 = 8192,
  MAIA_DEBUG_USER4 = 16384,
  MAIA_DEBUG_USER5 = 32768
};

//! The DEBUG class for the maia
/** \brief This class handles the DEBUG messaging for the maia.
 *
 * For debugging you mainly use two macros :
 * SET_DEBUG_LEVEL(debuglevel) and DEBUG( message, debuglevel)
 * The SET_DEBUG_LEVEL macro accepts arbitrary or combinations of
 * debug levels, so that you can display only the messages you want
 * to observe.
 * The message that is given to the DEBUG macro is only displayed, if
 * the connected debug level was given to the SET_DEBUG_LEVEL function
 * (as a single value or a combination).
 * If you temporarily want to suppress the DEBUG output you can use the
 * macros DEBUG_DISPLAY_ON and DEBUG_DISPLAY_OFF, which do what
 * their names imply.
 *
 * The common debug policy in the maia development is that you use
 * MAIA_DEBUG_TRACE only for marking the entry and return point of a
 * function. As long as your code is still experimental use the
 * MAIA_DEBUG_USER levels (choose the numbers at will). If you handover
 * your code, i.e. if it is added to the project all MAIA_DEBUG_USER
 * levels must be removed.
 * All DEBUG output that might be usefull later (i.e. also for other
 * programmers) must be set to MAIA_DEBUG_LEVEL levels.
 *
 * So (ideally speaking) you should never get code with MAIA_DEBUG_USER
 * levels set and you also should never hand over code which still
 * includes MAIA_DEBUG_USER levels.
 *
 * ---------------------------EXAMPLE :--------------------------------
 *
 * SET_DEBUG_LEVEL(MAIA_DEBUG_USER1|MAIA_DEBUG_USER3);
 *
 * DEBUG ( "Hello", MAIA_DEBUG_USER2 );
 * // this prints "Hello" if the DEBUG level is set to MAIA_DEBUG_USER2
 *
 * DEBUG ( "Solver::Solver" << " entry, solverId: " << m_solverId , MAIA_DEBUG_USER2 );
 * // hint: DEBUG works exactly like cout
 *
 * ----------------------- END OF EXAMPLE :----------------------------
 *
 *
 *  !!!!!!!!!!  Important note  !!!!!!!!!!!!
 *  SINCE DEBUG IS A MACRO, NEVER BREAK THE LINE IN WHICH YOU USE IT,
 *  I.E. DEBUG WON'T WORK WITH ARGUMENTS OF MULTIPLE LINES !
 */


class MDebug {
 public:
  static MString debugMessage;

  static void maiaerror(MString message, MDebugLevel debugLevel);
  static void maiaprint(const MString& message, MDebugLevel debugLevel);
  static void setLevel(MInt debugLevel);
  static void setLevelNotEnabled(const MInt debugLevel);
  static void displayOn();
  static void displayOff();

  static const MInt m_minLevel = MAIA_DEBUG_ASSERTION;
  static const MInt m_maxLevel = MAIA_DEBUG_USER5;

 private:
  static MInt m_debugLevel;
  static MBool m_debugOn;
  static std::basic_string<char> m_traceSpaces;
};

inline void MDebug::displayOn() { m_debugOn = true; }

inline void MDebug::displayOff() { m_debugOn = false; }

inline void MDebug::setLevel(MInt debugLevel) {
  displayOn();
  m_debugLevel = debugLevel;
}

inline void MDebug::setLevelNotEnabled(const MInt debugLevel) {
  if(debugLevel != 0) {
    cerr0 << "WARNING: trying to set debug level " << debugLevel << " but MAIA_DEBUG_FUNCTION is not enabled."
          << std::endl;
  }
}

inline void MDebug::maiaprint(const MString& message, MDebugLevel debugLevel) {
  // This function seperates the debug levels and
  // displays the messages which have the right debug level.
  if(m_debugOn) {
    auto a = (MUint)m_maxLevel;
    auto b = (MUint)m_debugLevel;
    ldiv_t adiv;
    do {
      adiv = ldiv(b, a);
      if(adiv.quot == 1 && a == (MUint)debugLevel) {
        if(debugLevel == MAIA_DEBUG_TRACE_OUT) {
          m_traceSpaces.erase(0, 2);
        }
        m_log << m_traceSpaces << "MDebug: "
              << "[" << globalDomainId() << "]:" << message << std::endl;
        if(debugLevel == MAIA_DEBUG_TRACE_IN) {
          m_traceSpaces += "  ";
        }
      }
      if(b >= a) {
        b -= a;
      }
      a = a / 2;
    } while(adiv.rem != 0); // repeat as long as there
    // is a rest and a level
  }
}

#if(defined(MAIA_DEBUG_FUNCTION) || defined(MAIA_PROFILING))
/// Helper class for automatic tracing
namespace maia {
namespace debug {
struct Tracer {
  Tracer(const MString& fun, [[maybe_unused]] const MString& loc, const MInt timingId, const MBool profile = false)
    : m_timingId(timingId), m_fun(fun), m_profile(profile) {
#ifdef MAIA_PROFILING
    if(profile) {
      PROFILE(timingId, fun, MAIA_DEBUG_TRACE_IN);
    } else {
      DEBUG(fun + " entry (" + loc + ")", MAIA_DEBUG_TRACE_IN);
    }
#else
    DEBUG(fun + " entry (" + loc + ")", MAIA_DEBUG_TRACE_IN);
#endif
  }
  ~Tracer() {
#ifdef MAIA_PROFILING
    if(m_profile) {
      PROFILE(m_timingId, m_fun, MAIA_DEBUG_TRACE_OUT);
    } else {
      DEBUG(m_fun << " return", MAIA_DEBUG_TRACE_OUT);
    }
#else
    DEBUG(m_fun << " return", MAIA_DEBUG_TRACE_OUT);
#endif
  }

 private:
  const MInt m_timingId = -1;
  const MString m_fun{};
  const MBool m_profile = false;
};
} // namespace debug
} // namespace maia
#endif

#endif // DEBUG_H
