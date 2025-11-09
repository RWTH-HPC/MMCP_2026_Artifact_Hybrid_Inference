// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TIMER_H
#define TIMER_H

#include <iomanip>
#include <string>
#include <time.h>
#include "COMM/mpioverride.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maialikwid.h"
#include "IO/infoout.h"
#include "compiler_config.h"
#include "functions.h"
#include "globalvariables.h"

#if defined(MAIA_MS_COMPILER)
#include <Windows.h>
#else
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#endif

// Note: required for writing to m_log
#ifndef PVPLUGIN
extern InfoOutFile m_log;
extern InfoOutFile maia_res;
#else
extern std::ostream& m_log;
extern std::ostream& maia_res;
#endif

#if defined(MAIA_MS_COMPILER)
static const unsigned __int64 epoch = ((unsigned __int64)116444736000000000ULL);
#endif

#ifdef MAIA_TIMER_FUNCTION
#define NEW_TIMER_GROUP(id, groupName) const MInt id = timers().newGroup(groupName)
#define NEW_TIMER(id, timerName, groupId) const MInt id = timers().newGroupTimer(timerName, groupId)
#define NEW_SUB_TIMER(id, timerName, timerId) const MInt id = timers().newSubTimer(timerName, timerId)
#define NEW_TIMER_GROUP_STATIC(id, groupName) static const MInt id = timers().newGroup(groupName)
#define NEW_TIMER_STATIC(id, timerName, groupId) static const MInt id = timers().newGroupTimer(timerName, groupId)
#define NEW_SUB_TIMER_STATIC(id, timerName, timerId) static const MInt id = timers().newSubTimer(timerName, timerId)
#define NEW_TIMER_GROUP_NOCREATE(id, groupName) id = timers().newGroup(groupName)
#define NEW_TIMER_NOCREATE(id, timerName, groupId) id = timers().newGroupTimer(timerName, groupId)
#define NEW_SUB_TIMER_NOCREATE(id, timerName, timerId) id = timers().newSubTimer(timerName, timerId)
#define RECORD_TIMER_START(timerId) timers().recordTimerStart(timerId, AT_)
#define RECORD_TIMER_STOP(timerId) timers().recordTimerStop(timerId, AT_)
#define RETURN_TIMER(timerId) timers().returnTimer(timerId, AT_)
#define RETURN_TIMER_TIME(timerId) timers().returnTimerTime(timerId)
#define STOP_ALL_TIMERS() timers().stopAllTimers(AT_)
#define RECORD_TIMER(timerId) timers().recordTimer(timerId, AT_)
#define RECORD_TIMERS() timers().recordTimers(AT_)
#define STOP_ALL_RECORD_TIMERS() timers().stopAllRecordTimers(AT_)
#define DISPLAY_TIMER(timerId) timers().displayTimer(timerId)
#define DISPLAY_TIMER_INTERM(timerId) timers().displayTimerNoToggleDisplayed(timerId)
#define DISPLAY_TIMER_OFFSET(timerId, ivl)                                                                             \
  if(globalTimeStep % ivl == 0) timers().displayTimerNoToggleDisplayed(timerId)
#define DISPLAY_ALL_GROUP_TIMERS(groupId) timers().displayAllTimers(groupId)
#define DISPLAY_ALL_TIMERS() timers().displayAllTimers()
#define RESET_TIMER(timerId) timers().resetTimer(timerId, AT_)
#define RESET_TIMERS() timers().resetTimers(AT_)
#define RESET_RECORD(timerId) timers().resetRecord(timerId)
#define SET_RECORD(timerId, timerValue) timers().resetRecord(timerId, timerValue)
#define RESET_ALL_RECORDS() timers().resetRecords()
#define SET_RECORD(timerId, timerValue) timers().resetRecord(timerId, timerValue)
#else
#define NEW_TIMER_GROUP(id, groupName)                                                                                 \
  do {                                                                                                                 \
  } while(false)
#define NEW_TIMER(id, timerName, groupId)                                                                              \
  do {                                                                                                                 \
  } while(false)
#define NEW_SUB_TIMER(id, timerName, timerId)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define NEW_TIMER_GROUP_STATIC(id, groupName)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define NEW_TIMER_STATIC(id, timerName, groupId)                                                                       \
  do {                                                                                                                 \
  } while(false)
#define NEW_SUB_TIMER_STATIC(id, timerName, timerId)                                                                   \
  do {                                                                                                                 \
  } while(false)
#define NEW_TIMER_GROUP_NOCREATE(id, groupName)                                                                        \
  do {                                                                                                                 \
  } while(false)
#define NEW_TIMER_NOCREATE(id, timerName, groupId)                                                                     \
  do {                                                                                                                 \
  } while(false)
#define NEW_SUB_TIMER_NOCREATE(id, timerName, timerId)                                                                 \
  do {                                                                                                                 \
  } while(false)
#define RECORD_TIMER_START(timerId)                                                                                    \
  do {                                                                                                                 \
  } while(false)
#define RECORD_TIMER_STOP(timerId)                                                                                     \
  do {                                                                                                                 \
  } while(false)
#define RETURN_TIMER(timerId)                                                                                          \
  do {                                                                                                                 \
  } while(false)
#define STOP_ALL_TIMERS()                                                                                              \
  do {                                                                                                                 \
  } while(false)
#define RECORD_TIMER(timerId)                                                                                          \
  do {                                                                                                                 \
  } while(false)
#define DISPLAY_TIMER(timerId)                                                                                         \
  do {                                                                                                                 \
  } while(false)
#define DISPLAY_ALL_GROUP_TIMERS(groupId)                                                                              \
  do {                                                                                                                 \
  } while(false)
#define DISPLAY_ALL_TIMERS()                                                                                           \
  do {                                                                                                                 \
  } while(false)
#define RESET_TIMER(timerId)                                                                                           \
  do {                                                                                                                 \
  } while(false)
#define RESET_TIMERS()                                                                                                 \
  do {                                                                                                                 \
  } while(false)
#define RESET_RECORD(timerId)                                                                                          \
  do {                                                                                                                 \
  } while(false)
#define RESET_ALL_RECORDS()                                                                                            \
  do {                                                                                                                 \
  } while(false)
#define SET_RECORD(timeValue, timerId)                                                                                 \
  do {                                                                                                                 \
  } while(false)
#endif

/// \brief MTimers manages all MAIA Timers and allows primitive profiling.
///
/// Usage:
/// - NEW_TIMER(string name) creates a new timer with name "name" and
/// returns its index (a MInt that you can use to access it).
/// - RESET_TIMER(MInt timerId): resets timerId.
/// - START_TIMER(timerId)/STOP_TIMER(timerId) work as expected.
/// - DISPLAY_TIMER(timerId) writes the timerId name and time to the log.
/// - DISPLAY_ALL_TIMERS: displays all timers.
///
/// Example:
///
/// RESET_TIMER;
/// MInt myTimer;
/// myTimer = NEW_TIMER("My timer");
///
/// START_TIMER(myTimer);
/// f1(); // Function will be timed.
/// STOP_TIMER(myTimer);
/// f2(); // Function will not be timed.
/// START_TIMER(myTimer);
/// f3(); // Function will be timed.
/// STOP_TIMER(myTimer);
///
/// DISPLAY_TIMER(myTimer);
///
/// NOTE: MAIA_TIMER_FUNCTION macro enable/disables timing MAIA.
///
class MTimers {
  friend MTimers& timers();

 public:
  inline MInt newGroup(const std::string groupName);
  inline MInt newGroupTimer(const std::string timerName, const MInt groupId);
  inline MInt newSubTimer(const std::string timerName, const MInt timerId);
  inline MFloat returnTimer(const MInt timerId, const MString pos);
  inline MFloat returnTimerTime(const MInt timerId);
  inline void resetTimer(const MInt timerId, const MString pos);       // reset
  inline void recordTimer(const MInt timerId, const MString pos);      // record
  inline void recordTimerStart(const MInt timerId, const MString pos); // reset + start
  inline void recordTimerStop(const MInt timerId, const MString pos);  // stop + record
  inline void recordTimers(const MString pos);
  inline void stopAllTimers(const MString pos);
  inline void stopAllRecordTimers(const MString pos); // stop all + record stopped
  inline void displayTimer(const MInt timerId);
  inline void displayTimerNoToggleDisplayed(const MInt timerId);
  inline void displayAllTimers();
  inline void displayAllTimers(const MInt groupId);
  inline void resetTimers(const MString pos);
  inline void resetRecord(const MInt timerId, const MFloat timerValue = 0.0);
  inline void resetRecords();
  inline MInt isRunning(const MInt timerId) const { return m_timers[timerId].status == Timer::Running; };

 private:
  inline void startTimer(const MInt timerId, const MString pos);
  inline void stopTimer(const MInt timerId, const MString pos);

 private:
  MTimers() {
#if defined(MAIA_TIMER_FUNCTION)
    m_log << "MTimers: timers enabled." << std::endl;
#endif
#if defined(MAIA_TIMER_DEBUG)
    m_log << "MTimers: timer debug enabled." << std::endl;
#endif
#if defined(MAIA_TIMER_CHECKS)
    m_log << "MTimers: timer checks enabled." << std::endl;
#endif
#if defined(MAIA_ASSERT_TIMER_CHECKS)
    m_log << "MTimers: assert timer checks enabled." << std::endl;
#endif
  }
  ~MTimers() {}
  // delete: copy construction, and copy assignment
  MTimers(MTimers&);
  MTimers& operator=(const MTimers&);

  struct Timer {
    Timer(const MString n, const MInt g, const MInt id, const MInt p)
      : name(n),
        group(g),
        timerId(id),
        parent(p),
        cpuTime(0),
        oldCpuTime(0),
        recordedTime(0),
        status(Timer::Uninitialized),
        subTimers(0) {}
    MString name;    ///< Timer Name
    MInt group = -1; ///< Group Id
    MInt timerId = -1;
    MInt parent = -1;    ///< Parent timer id
    MFloat cpuTime = 0.0;      ///< CPU time
    MFloat oldCpuTime = 0.0;   ///< Old CPU time (for timer restart)
    MFloat recordedTime = 0.0; ///< Time recorded on the timer.
    MInt status = -1;          ///< Timer's status, see enum:
    enum { Uninitialized = 0, Running = 1, Stopped = 2 };
    std::vector<MInt> subTimers{};
    MBool displayed = false;
  };

  std::vector<std::string> m_groups;
  std::vector<Timer> m_timers;

  inline MFloat time();
  inline void displayTimer_(const MInt timerId,
                            const MBool toggleDisplayed = true,
                            const MInt tIndent = 0,
                            const MFloat superTime = -1.0);
  inline void displayTimerHeader_();
  inline void displayTimerGroupHeader_(const MInt groupId);
  inline MInt indent(const MInt pIndent) const { return pIndent + 2; };
};

MTimers& timers();

/// Returns Wall-Clock time in seconds
inline MFloat MTimers::time() {
  /// Timers are not portable.
  /// We should use the C++ <chrono> library instead.
  /// The factors convert the cpu time to seconds
#ifdef MAIA_MPI_TIMER
  return MPI_Wtime();
#else
#if _POSIX_TIMERS > 0
  timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return static_cast<MFloat>(t.tv_sec) + static_cast<MFloat>(t.tv_nsec / 1000000000.0);
#else
  struct timeval t;
  gettimeofday(&t, nullptr);
  return static_cast<MFloat>(t.tv_sec) + static_cast<MFloat>(t.tv_usec / 1000000.0);
#endif
#endif
}

inline void MTimers::resetRecord(const MInt timerId, const MFloat timerValue) {
  m_timers[timerId].recordedTime = timerValue;
}

inline void MTimers::resetRecords() {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    resetRecord((MInt)timerId);
  }
}

inline void MTimers::recordTimerStop(const MInt timerId, const MString pos) {
  if(timerId < 0) {
    return;
  }
#ifdef MAIA_SYNCHRONIZE_TIMERS
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  stopTimer(timerId, pos);
  recordTimer(timerId, pos);
}

inline void MTimers::recordTimer(const MInt timerId, const MString pos) {
  m_timers[timerId].recordedTime += returnTimer(timerId, pos);
}

inline void MTimers::recordTimers(const MString pos) {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    recordTimer(timerId, pos);
  }
}

inline void MTimers::resetTimer(const MInt timerId, [[maybe_unused]] const MString pos) {
#if defined(MAIA_TIMER_CHECKS)
  if(m_timers[timerId].status == Timer::Running) {
    const MString msg = "The timer #" + std::to_string(timerId) + " '" + m_timers[timerId].name
                        + "' can't be reset because it is running! " + pos;
    m_log << msg << std::endl;
    std::cerr << msg << std::endl;
#if defined(MAIA_ASSERT_TIMER_CHECKS)
    TERMM(1, msg);
#endif
  }
#endif

  m_timers[timerId].cpuTime = 0.0;
  m_timers[timerId].status = Timer::Stopped;
}

inline void MTimers::resetTimers(const MString pos) {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    resetTimer(timerId, pos);
  }
}

/// Creates a new timer group and returns its groupId.
inline MInt MTimers::newGroup(const std::string name) {
  m_groups.push_back(name);
#if defined(MAIA_TIMER_DEBUG)
  m_log << "New timer group " << m_groups.size() - 1 << " '" << name << "' " << std::endl;
#endif
  return m_groups.size() - 1;
}

/// Creates a new timer and returns its timerId.
inline MInt MTimers::newGroupTimer(const std::string name, const MInt groupId) {
  ASSERT(static_cast<std::size_t>(groupId) < m_groups.size() && groupId > -1,
         "groupId: " << groupId << " does not exists | name: " << name);
  const MInt newTimerId = m_timers.size();
  m_timers.push_back(Timer(name, groupId, newTimerId, -1));
#if defined(MAIA_TIMER_DEBUG)
  m_log << "New group timer " << newTimerId << " '" << name << "' " << groupId << std::endl;
#endif
  return newTimerId;
}

/// Creates a new timer and returns its timerId.
inline MInt MTimers::newSubTimer(const std::string name, const MInt timerId) {
  if(timerId < 0) {
#if defined(MAIA_TIMER_DEBUG)
    m_log << "cannot create subTimer '" << name << "'" << std::endl;
#endif
    return -1;
  }

  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(),
         "timerId " << timerId << " does not exist when trying to create subtimer with name " << name);

  const MInt groupId = m_timers[timerId].group;
  const MInt newTimerId = m_timers.size();
  m_timers.push_back(Timer(name, groupId, newTimerId, timerId));
#if defined(MAIA_TIMER_DEBUG)
  m_log << "New subtimer " << newTimerId << " '" << name << "' " << groupId << " " << timerId << std::endl;
#endif
  m_timers[timerId].subTimers.push_back(newTimerId);
  return newTimerId;
}

inline void MTimers::recordTimerStart(const MInt timerId, const MString pos) {
  if(timerId < 0) {
    return;
  }
#ifdef MAIA_SYNCHRONIZE_TIMERS
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  resetTimer(timerId, pos);
  startTimer(timerId, pos);
}

inline void MTimers::startTimer(const MInt timerId, [[maybe_unused]] const MString pos) {
  ASSERT(m_timers[timerId].status != Timer::Running,
         "The timer " << m_timers[timerId].name << " with id: " << timerId
                      << " can't be started because it is already running! " << pos);

#if defined(MAIA_TIMER_CHECKS)
  const MInt parent = m_timers[timerId].parent;
  if(parent != -1) {
    if(m_timers[parent].status != Timer::Running) {
      const MString msg = "The timer #" + std::to_string(timerId) + " '" + m_timers[timerId].name
                          + "' can't be started because its parent timer #" + std::to_string(parent) + " '"
                          + m_timers[parent].name + "' is not running! " + pos;
      m_log << msg << std::endl;
      std::cerr << msg << std::endl;
#if defined(MAIA_ASSERT_TIMER_CHECKS)
      TERMM(1, msg);
#endif
    }
  }
#endif

#if defined(MAIA_TIMER_DEBUG)
  m_log << "start timer #" + std::to_string(timerId) + " '" + m_timers[timerId].name + "' " << pos << std::endl;
#endif

  const MFloat t = time();
  m_timers[timerId].oldCpuTime = m_timers[timerId].cpuTime;
  m_timers[timerId].cpuTime = t;
  m_timers[timerId].status = Timer::Running;

  // Enable likwid counter if likwid is enabled
#ifdef WITH_LIKWID
  LIKWID_MARKER_START(std::to_string(timerId).c_str());
#endif
}


/// Returns the timer Value.
inline MFloat MTimers::returnTimer(const MInt timerId, [[maybe_unused]] const MString pos) {
  if(timerId < 0) {
    TERMM(1, "Invalid timer id");
    return 0.0;
  }

#if defined(MAIA_TIMER_CHECKS)
  if(m_timers[timerId].status != Timer::Stopped) {
    const MString msg = "The timer #" + std::to_string(timerId) + " '" + m_timers[timerId].name
                        + "' needs to be stopped before returnTimer() is called! " + pos;
    m_log << msg << std::endl;
    std::cerr << msg << std::endl;
#if defined(MAIA_ASSERT_TIMER_CHECKS)
    TERMM(1, msg);
#endif
  }
#endif

#if MAIA_TIMERS_AVERAGE_OVER_DOMAINS
  const MFloat t = m_timers[timerId].cpuTime;
  MFloat tmp_rcv = 0.0;
  MPI_Reduce(&t, &tmp_rcv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  return tmp_rcv / globalNoDomains();
#else
  return m_timers[timerId].cpuTime;
#endif
}

// Returns the recorded time
inline MFloat MTimers::returnTimerTime(const MInt timerId) { return m_timers[timerId].recordedTime; }

/// Stops the timer and sets its final value.
inline void MTimers::stopTimer(const MInt timerId, [[maybe_unused]] const MString pos) {
  if(timerId < 0) {
#if defined(MAIA_TIMER_DEBUG)
    m_log << "cannot stop timer " << timerId << " " << pos << std::endl;
#endif
    return;
  }

  if(m_timers[timerId].status == Timer::Running) {
    const MFloat t = time();
    m_timers[timerId].cpuTime = t - m_timers[timerId].cpuTime + m_timers[timerId].oldCpuTime;
    m_timers[timerId].status = Timer::Stopped;

#if defined(MAIA_TIMER_CHECKS)
    for(auto subTimerId : m_timers[timerId].subTimers) {
      if(m_timers[subTimerId].status == Timer::Running) {
        const MString msg = "The timer #" + std::to_string(timerId) + " '" + m_timers[timerId].name
                            + "' can't be stopped because its sub-timer #" + std::to_string(subTimerId) + " '"
                            + m_timers[subTimerId].name + "' is still running! " + pos;
        m_log << msg << std::endl;
        std::cerr << msg << std::endl;
#if defined(MAIA_ASSERT_TIMER_CHECKS)
        TERMM(1, msg);
#endif
      }
    }
#endif

#if defined(MAIA_TIMER_DEBUG)
    m_log << "stop timer #" + std::to_string(timerId) + " '" + m_timers[timerId].name + "' " << pos << std::endl;
#endif

    // Stop likwid counter if likwid is enabled
#ifdef WITH_LIKWID
    LIKWID_MARKER_STOP(std::to_string(timerId).c_str());
#endif
  } else {
#if defined(MAIA_TIMER_CHECKS)
    const MString msg = "The timer '" + m_timers[timerId].name + "' can't be stopped because it is not running! " + pos;
    m_log << msg << std::endl;
    std::cerr << msg << std::endl;
#if defined(MAIA_ASSERT_TIMER_CHECKS)
    TERMM(1, msg);
#endif
#endif
  }
}


// Stops all timers.
inline void MTimers::stopAllTimers(const MString pos) {
  for(std::size_t i = 0, e = m_timers.size(); i != e; ++i) {
    if(m_timers[i].status == Timer::Running) {
      stopTimer(i, pos);
    }
  }
}

// Stops all timers and record the timers that were stopped
inline void MTimers::stopAllRecordTimers(const MString pos) {
  // for(std::size_t i = 0, e = m_timers.size(); i != e; ++i) {
  for(MInt i = m_timers.size() - 1; i >= 0; i--) {
    if(m_timers[i].status == Timer::Running) {
#ifdef MAIA_SYNCHRONIZE_TIMERS
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      stopTimer(i, pos);
      recordTimer(i, pos);
    }
  }
}


inline void MTimers::displayTimer_(const MInt timerId, const MBool toggleDisplayed, const MInt tIndent,
                                   const MFloat superTime) {
  if(m_timers[timerId].displayed) {
    return;
  }

  /* MBool running = false; */
  if(m_timers[timerId].status == Timer::Running) {
    // NOTE: test this if needed, but not used at the moment
    TERMM(1, "Timer is still running");
    /* running = true; */
    /* stopTimer(timerId, AT_); */
  }
  m_log.width(50);
  m_log.setf(std::ios::left);
  std::stringstream indentedName;

  // Calculate time relative to the parent timer
  MFloat percentage;
  if(superTime < F0) {
    // If the parent time is less than zero, that means that there is no parent timer
    // and the percentage should be 100%
    percentage = 100.0;
  } else if(approx(superTime, F0, MFloatEps)) {
    // If the parent time is approximately zero, that probably means that the timer was never
    // run - therefore the percentage is set to 0%
    percentage = F0;
  } else {
    // Otherwise calculate the percentage as the fraction of this timer vs. the parent timer times 100%
    percentage = 100.0 * m_timers[timerId].recordedTime / superTime;
  }

  indentedName << std::string(tIndent, ' ');
  indentedName << "[" << std::fixed << std::setprecision(1) << std::setw(4) << std::setfill('0') << std::right
               << percentage << std::left << "%] ";
  indentedName << m_timers[timerId].name;
  m_log << indentedName.str() << std::right;
  m_log.precision(6);
  m_log.width(20);
  m_log << m_timers[timerId].recordedTime << std::left << " [sec]";
  // Show output of likwid performance counters if likwid is enabled
#ifdef WITH_LIKWID
  // If the timer wasn't called set the MFLops to 0.00
  m_log << "   ";
  if(approx(m_timers[timerId].recordedTime, F0, MFloatEps)) {
    m_log << "0.00";
  } else {
    m_log << "${timer_" << timerId << "}";
  }
  m_log << " [DP MFlops/s]";
#endif
  if(toggleDisplayed) m_timers[timerId].displayed = true;
  m_log << std::endl;
  for(std::size_t sub = 0, last = m_timers[timerId].subTimers.size(); sub < last; ++sub) {
    const MInt new_indent = indent(tIndent);
    displayTimer_(m_timers[timerId].subTimers[sub], toggleDisplayed, new_indent, m_timers[timerId].recordedTime);
  }
  /* if(running) { */
  /*   startTimer(timerId, AT_); */
  /* } */
}

inline void MTimers::displayTimerHeader_() {}

inline void MTimers::displayTimerGroupHeader_(const MInt groupId) {
  m_log << "--------------------------------------------------------------------------------" << std::endl;
  m_log.width(50);
  m_log.precision(12);
  m_log.setf(std::ios::left);
  m_log << "Group";
  m_log.width(40);
  m_log << m_groups[groupId] << std::endl;
}

inline void MTimers::displayAllTimers() {
  ASSERT(m_timers.size() > 0, "ERROR: no timers have been created!");
  for(std::size_t groupId = 0, e = m_groups.size(); groupId != e; ++groupId) {
    displayAllTimers(groupId);
  }
}

inline void MTimers::displayAllTimers(const MInt groupId) {
  ASSERT(m_timers.size() > 0, "ERROR: no timers have been created!");
  ASSERT(static_cast<std::size_t>(groupId) < m_groups.size() && groupId > -1, "ERROR: groupId does not exists");
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    m_timers[timerId].displayed = false;
  }
  displayTimerGroupHeader_(groupId);
  displayTimerHeader_();
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    if(m_timers[timerId].group == groupId) {
      displayTimer_(timerId);
    }
  }
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    m_timers[timerId].displayed = false;
  }
}

inline void MTimers::displayTimer(const MInt timerId) {
  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(), "ERROR: timer timerId does not exist");
  displayTimerHeader_();
  displayTimer_(timerId);
}

inline void MTimers::displayTimerNoToggleDisplayed(const MInt timerId) {
  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(), "ERROR: timer timerId does not exist");
  displayTimerHeader_();
  displayTimer_(timerId, false);
}

//------------------------------------------------------------------------------

class Profile;
class FunctionTiming;

/**
 * \brief This class counts the static execution time of a function
 * \author Lennart Schneiders
 * \date 14.02.2013
 */
class FunctionTiming {
 public:
  explicit FunctionTiming(std::string name);
  ~FunctionTiming();
  FunctionTiming& operator=(const FunctionTiming& t);
  void in();
  void out();
  MFloat getInitCpuTime() const { return m_initCpuTime; }
  MFloat getDeltaCpuTime() const { return m_deltaCpuTime; }
  MFloat getInitWallTime() const { return m_initWallTime; }
  MFloat getDeltaWallTime() const { return m_deltaWallTime; }
  MString getName() const { return m_name; }

 private:
  MFloat m_initCpuTime = -1.0;
  MFloat m_deltaCpuTime = -1.0;
  MFloat m_tmpCpuTime = -1.0;
  MFloat m_initWallTime = -1.0;
  MFloat m_deltaWallTime = -1.0;
  MFloat m_tmpWallTime = -1.0;
  MString m_name = "";
};

MBool operator<(const FunctionTiming& a, const FunctionTiming& b);

/**
 * \brief This class collects all function timings and produces a profiling for certain areas of the code
 * \author Lennart Schneiders
 * \date 14.02.2013
 */
class Profile {
 public:
  explicit Profile(const std::string& name) : m_initCpuTime(cpuTime()), m_initWallTime(wallTime()), m_name(name) {}
  ~Profile();
  static MInt getTimingId(std::string name);
  static MString printTime(MFloat secs);
  static std::vector<FunctionTiming> s_functionTimings;

 private:
  const MFloat m_initCpuTime;
  const MFloat m_initWallTime;
  const MString m_name;
};


void logDuration_(const MFloat timeStart,
                  const MString module,
                  const MString comment,
                  const MPI_Comm comm,
                  const MInt domainId,
                  const MInt noDomains);

void logDurations(std::vector<std::pair<MFloat, MString>>& durations, const MString module, const MPI_Comm comm,
                  const MInt domainId, const MInt noDomains);

#endif // TIMER_H
