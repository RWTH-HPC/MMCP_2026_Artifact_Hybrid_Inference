// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DLBTIMER_H
#define DLBTIMER_H

#include <iomanip>
#include <sys/times.h>
#include <unistd.h>
#include "INCLUDE/maiatypes.h"
#include "IO/infoout.h"
#include "UTIL/timer.h"
#include "functions.h"
#include "globalvariables.h"

// Note: required for writing to m_log, since dlbtimer.h is itself included in
// globalvariables.h
#ifndef PVPLUGIN
extern InfoOutFile m_log;
extern InfoOutFile maia_res;
#else
extern std::ostream& m_log;
extern std::ostream& maia_res;
#endif

namespace maia {
namespace dlbTimer {

class DlbTimer {
 private:
  /// Unique DLB timer id among all DLB timers
  const MInt m_dlbTimerId = -1;

  /// Timer indices
  const MInt m_loadTimerId = 0;
  const MInt m_idleTimerId = 1;

  /// Stores if the timers are enabled
  MBool m_dlbTimersEnabled = false;
  /// Current time when the timers were started
  MFloat m_timerWallTime[2] = {0.0, 0.0}; // wall time
  MFloat m_timerUserTime[2] = {0.0, 0.0}; // process user time
  /// Recorded time of the timers
  MFloat m_recordedWallTime[2] = {0.0, 0.0};
  MFloat m_recordedUserTime[2] = {0.0, 0.0};
  /// Status of the timers
  MBool m_timerRunning[2] = {false, false};

 public:
  /// Constructor, just set the timer id
  DlbTimer(const MInt timerId) : m_dlbTimerId(timerId) {}

  /// Enable timers for dynamic load balancing, disabled by default since not supported by all
  /// solvers/run loops
  inline void enableDlbTimers() {
    if(m_dlbTimersEnabled) {
      TERMM(1, "Error: load balancing timers already enabled.");
    }
    m_dlbTimersEnabled = true;

    // Reset timer status
    for(MBool& i : m_timerRunning) {
      i = false;
    }
  }

  /// Temporarily disable timers, e.g. during adaptation
  inline void disableDlbTimers() {
    if(!m_dlbTimersEnabled) {
      TERMM(1, "Error: load balancing timers not enabled.");
    }
    m_dlbTimersEnabled = false;
  }

  inline void reEnableDlbTimer() {
    if(m_dlbTimersEnabled) {
      TERMM(1, "Error: load balancing timers already enabled.");
    }
    m_dlbTimersEnabled = true;
  }

  /// \brief Return if timers are enabled
  inline MBool dlbTimersEnabled() const { return m_dlbTimersEnabled; }

  inline MBool isLoadTimerRunning() const { return m_timerRunning[m_loadTimerId]; }

  /// \brief Start the load timer
  inline void startLoadTimer(const MString& name) {
    if(!m_dlbTimersEnabled) {
      return;
    }
    if(m_timerRunning[m_idleTimerId]) {
      TERMM(1, name + "; Error: cannot start load timer while idle timer still running.");
    }
    startSolverTimer(m_loadTimerId, name);
  }

  /// \brief Stop the load timer
  inline void stopLoadTimer(const MString& name) {
    if(!m_dlbTimersEnabled) {
      return;
    }
    stopSolverTimer(m_loadTimerId, name);
  }

  /// \brief Reset the load record
  inline void resetLoadRecord() { resetSolverTimer(m_loadTimerId); }

  /// \brief Return the load record
  inline MFloat returnLoadRecord(const MInt mode = 0) const {
#ifdef MAIA_DEBUG_DLB_TIMER
    const MFloat diff = m_recordedWallTime[m_loadTimerId] - m_recordedUserTime[m_loadTimerId];
    m_log << "DLB timer #" << m_dlbTimerId << " "
          << " returnLoadRecord: " << m_recordedWallTime[m_loadTimerId] << " user " << m_recordedUserTime[m_loadTimerId]
          << " diff " << diff << std::endl;
#endif
    return (mode == 0) ? m_recordedWallTime[m_loadTimerId] : m_recordedUserTime[m_loadTimerId];
  }

  /// \brief Start the idle timer
  inline void startIdleTimer(const MString& name) {
    if(!m_dlbTimersEnabled) {
      return;
    }
    if(m_timerRunning[m_loadTimerId]) {
      TERMM(1, name + "; Error: cannot start idle timer while load timer still running.");
    }
    startSolverTimer(m_idleTimerId, name);
  }

  /// \brief Stop the idle timer
  inline void stopIdleTimer(const MString& name) {
    if(!m_dlbTimersEnabled) {
      return;
    }
    stopSolverTimer(m_idleTimerId, name);
  }

  /// \brief Reset the idle record
  inline void resetIdleRecord() { resetSolverTimer(m_idleTimerId); }

  /// \brief Return the idle record
  inline MFloat returnIdleRecord(const MInt mode = 0) const {
#ifdef MAIA_DEBUG_DLB_TIMER
    const MFloat diff = m_recordedWallTime[m_idleTimerId] - m_recordedUserTime[m_idleTimerId];
    m_log << "DLB timer #" << m_dlbTimerId << " "
          << " returnIdleRecord: " << m_recordedWallTime[m_idleTimerId] << " user " << m_recordedUserTime[m_idleTimerId]
          << " diff " << diff << std::endl;
#endif
    return (mode == 0) ? m_recordedWallTime[m_idleTimerId] : m_recordedUserTime[m_idleTimerId];
  }

  static MInt noTimers() {
    return 2; // Default: load and idle timer
  }

 private:
  // General DLB timer methods

  /// \brief Start the timer with the given id
  void startSolverTimer(const MInt timerId, const MString& name) {
    if(m_timerRunning[timerId]) {
      TERMM(1, name + "; error in startSolverTimer(" + std::to_string(timerId) + "): timer already running, dlbTimerId "
                   + std::to_string(m_dlbTimerId));
    } else {
      m_timerWallTime[timerId] = MPI_Wtime();
      m_timerUserTime[timerId] = cpuTime();

      m_timerRunning[timerId] = true;
#ifdef MAIA_DEBUG_DLB_TIMER
      m_log << "DLB timer #" << m_dlbTimerId << " startSolverTimer #" << timerId << " at time "
            << m_timerWallTime[timerId] << " from " << name << std::endl;
#endif
    }
  }

  /// \brief Stop the timer with the given id
  void stopSolverTimer(const MInt timerId, const MString& name) {
    if(!m_timerRunning[timerId]) {
      TERMM(1, name + "; error in stopSolverTimer(" + std::to_string(timerId) + "): timer not running, dlbTimerId "
                   + std::to_string(m_dlbTimerId) + " DLB status is " + std::to_string(dlbTimersEnabled()));
    } else {
      const MFloat t_stop = MPI_Wtime();
      const MFloat t_diff = t_stop - m_timerWallTime[timerId];
      m_recordedWallTime[timerId] += t_diff;

      const MFloat t_user = cpuTime();
      const MFloat t_user_diff = t_user - m_timerUserTime[timerId];
      m_recordedUserTime[timerId] += t_user_diff;

      m_timerRunning[timerId] = false;
#ifdef MAIA_DEBUG_DLB_TIMER
      m_log << "DLB timer #" << m_dlbTimerId << " stopSolverTimer #" << timerId << " at time "
            << m_timerWallTime[timerId] << " with t_diff " << t_diff << " user time diff " << t_user_diff << " from "
            << name << std::endl;
#endif
    }
  }

  /// \brief Reset the timer with the given id
  void resetSolverTimer(const MInt timerId) {
    if(m_timerRunning[timerId]) {
      TERMM(1, "error in resetSolverTimer(" + std::to_string(timerId) + "): timer still running, dlbTimerId "
                   + std::to_string(m_dlbTimerId));
    } else {
      m_recordedWallTime[timerId] = 0.0;
      m_recordedUserTime[timerId] = 0.0;
#ifdef MAIA_DEBUG_DLB_TIMER
      m_log << "DLB timer #" << m_dlbTimerId << " resetSolverTimer #" << timerId << std::endl;
#endif
    }
  }
};

} // namespace dlbTimer
} // namespace maia


/// \brief Controller class for all DLB timers
class DlbTimerController {
#define ENSURE_VALID_TIMERID(timerId, at)                                                                              \
  do {                                                                                                                 \
    ASSERT(timerId >= 0 && timerId < noDlbTimers(), "invalid dlbTimerId: " + std::to_string(timerId));                 \
  } while(false)
#define ENSURE_VALID_TIMERMODE(mode, at)                                                                               \
  do {                                                                                                                 \
    ASSERT(mode >= 0 && mode <= 1, "invalid timer mode: " + std::to_string(mode));                                     \
  } while(false)

 public:
  /// \brief Create the given number of DLB timers
  void createDlbTimers(const MInt noTimers, const MBool ignore = false) {
    if(noDlbTimers() > 0) {
      TERMM(1, "createDlbTimers should only be called once!");
    }

    // High-resolution per-process CPU timer
    timespec tp;
    clock_getres(CLOCK_PROCESS_CPUTIME_ID, &tp);
    m_log << "DLB timer: CLOCK PROCESS CPUTIME RESOLUTION " << tp.tv_sec << "s " << tp.tv_nsec << "nsec" << std::endl;

    for(MInt i = 0; i < noTimers; i++) {
      const MInt newTimerId = m_dlbTimers.size();
      m_dlbTimers.emplace_back(newTimerId);
    }

    m_ignoreDlbTimers = ignore;
  }

  /// \brief Enable the given DLB timer
  void enableDlbTimers(const MInt dlbTimerId) {
    if(m_ignoreDlbTimers) {
      return;
    }

    ENSURE_VALID_TIMERID(dlbTimerId, AT_);
    m_dlbTimers[dlbTimerId].enableDlbTimers();
    m_enabled = true; // At least this one timer is enabled
  }

  /// \brief Enable all DLB timers (or those given by the array wasEnabled)
  void enableAllDlbTimers(const MBool* const wasEnabled = nullptr) {
    if(m_ignoreDlbTimers) {
      return;
    }

    for(MInt i = 0; i < noDlbTimers(); i++) {
      // Check if a status array is given and restore previous state
      MBool enableTimers = true;
      if(wasEnabled != nullptr) {
        enableTimers = wasEnabled[i];
      }

      if(enableTimers) {
        enableDlbTimers(i);
      }
    }
  }

  /// \brief Disable the given DLB timer
  void disableDlbTimers(const MInt dlbTimerId) {
    if(m_ignoreDlbTimers) {
      return;
    }

    ENSURE_VALID_TIMERID(dlbTimerId, AT_);
    m_dlbTimers[dlbTimerId].disableDlbTimers();

    MBool anyEnabled = false;
    for(MInt i = 0; i < noDlbTimers(); i++) {
      anyEnabled |= m_dlbTimers[dlbTimerId].dlbTimersEnabled();
    }
    m_enabled = anyEnabled;
  }

  /// \brief Disable all (enabled) DLB timers
  void disableAllDlbTimers(MBool* const wasEnabled = nullptr) {
    if(m_ignoreDlbTimers) {
      return;
    }

    if(m_runningTimerId != -1) {
      TERMM(1, "Cannot disable all DLB timers, timer " + std::to_string(m_runningTimerId) + " still running.");
    }
    for(MInt i = 0; i < noDlbTimers(); i++) {
      const MBool timersEnabled = dlbTimersEnabled(i);

      // Store current timer status if requested (pointer != nullptr passed to function)
      if(wasEnabled != nullptr) {
        wasEnabled[i] = timersEnabled;
      }

      if(timersEnabled) {
        disableDlbTimers(i);
      }
    }
    m_enabled = false; // All timers are disabled
  }

  /// \brief Return if the given DLB timer is enabled
  MBool dlbTimersEnabled(const MInt dlbTimerId) {
    if(m_ignoreDlbTimers) {
      return false;
    }

    ENSURE_VALID_TIMERID(dlbTimerId, AT_);
    return m_dlbTimers[dlbTimerId].dlbTimersEnabled();
  }

  /// \brief Start the load timer for the given DLB timer id
  void startLoadTimer(const MInt dlbTimerId, const MString& name) {
    if(m_ignoreDlbTimers) {
      return;
    }
    ENSURE_VALID_TIMERID(dlbTimerId, name);

    if(!dlbTimersEnabled(dlbTimerId)) {
      return;
    }

    // Check that no load/idle timer is currently running
    if(m_runningTimerId != -1) {
      TERMM(1, "Cannot start load timer " + std::to_string(dlbTimerId) + ", timer " + std::to_string(m_runningTimerId)
                   + " already running.");
    }

    m_dlbTimers[dlbTimerId].startLoadTimer(name);
    m_runningTimerId = dlbTimerId;
  }

  /// \brief Stop the load timer for the given DLB timer id
  void stopLoadTimer(const MInt dlbTimerId, const MString& name) {
    if(m_ignoreDlbTimers) {
      return;
    }
    ENSURE_VALID_TIMERID(dlbTimerId, name);

    if(!dlbTimersEnabled(dlbTimerId)) {
      return;
    }

    ASSERT(dlbTimerId == m_runningTimerId, "timer id does not match the running timer id");
    m_dlbTimers[dlbTimerId].stopLoadTimer(name);
    m_runningTimerId = -1;
  }

  /// \brief Start the idle timer for the given DLB timer id
  void startIdleTimer(const MInt dlbTimerId, const MString& name) {
    if(m_ignoreDlbTimers) {
      return;
    }
    ENSURE_VALID_TIMERID(dlbTimerId, name);

    if(!dlbTimersEnabled(dlbTimerId)) {
      return;
    }

    // Check that no load/idle timer is currently running
    if(m_runningTimerId != -1) {
      TERMM(1, "Cannot start idle timer " + std::to_string(dlbTimerId) + ", timer " + std::to_string(m_runningTimerId)
                   + " already running.");
    }

    m_dlbTimers[dlbTimerId].startIdleTimer(name);
    m_runningTimerId = dlbTimerId;
  }

  /// \brief Stop the idle timer for the given DLB timer id
  void stopIdleTimer(const MInt dlbTimerId, const MString& name) {
    if(m_ignoreDlbTimers) {
      return;
    }
    ENSURE_VALID_TIMERID(dlbTimerId, name);

    if(!dlbTimersEnabled(dlbTimerId)) {
      return;
    }

    ASSERT(dlbTimerId == m_runningTimerId, "timer id does not match the running timer id");
    m_dlbTimers[dlbTimerId].stopIdleTimer(name);
    m_runningTimerId = -1;
  }

  /// \brief Stop the currently running load timer and start the corresponding idle timer
  void stopLoadStartIdleTimer(const MString& name) {
    if(!m_enabled) {
      return; // If no timer is enabled there is nothing to do
    }

    if(m_runningTimerId < 0) {
      TERMM(1, "Cannot stop load and start idle timer, the running timer is: " + std::to_string(m_runningTimerId) + '/'
                   + name);
    }
    const MInt timerId = m_runningTimerId; // m_runningTimerId is reset in stopLoadTimer()
    stopLoadTimer(timerId, name);
    startIdleTimer(timerId, name);
  }

  /// \brief Stop the currently running idle timer and start the corresponding load timer
  void stopIdleStartLoadTimer(const MString& name) {
    if(!m_enabled) {
      return; // If no timer is enabled there is nothing to do
    }

    if(m_runningTimerId < 0) {
      TERMM(1, "Cannot stop idle and start load timer, the running timer is: " + std::to_string(m_runningTimerId));
    }
    const MInt timerId = m_runningTimerId; // m_runningTimerId is reset in stopIdleTimer()
    stopIdleTimer(timerId, name);
    startLoadTimer(timerId, name);
  }

  /// \brief Return if a timer is running
  MBool isTimerRunning() const { return (m_runningTimerId != -1); }

  MInt whichTimerIsRunning() const { return m_runningTimerId; }

  MBool isLoadTimerRunning(const MInt dlbTimerId) {
    return m_dlbTimers[dlbTimerId].isLoadTimerRunning() && m_runningTimerId == dlbTimerId;
  }

  void reEnableDlbTimer(const MInt dlbTimerId) { m_dlbTimers[dlbTimerId].reEnableDlbTimer(); }

  /// \brief Check the timer status during IO (no timer running and timers not enabled)
  void checkIOTimerStatus(const MString& name) const {
    TERMM_IF_COND(m_enabled || isTimerRunning(), "Timers still enabled and/or timer running " + name);
  }

  /// \brief Reset the records of all DLB timers
  void resetRecords() {
    if(m_ignoreDlbTimers) {
      return;
    }
    for(MInt i = 0; i < noDlbTimers(); i++) {
      m_dlbTimers[i].resetLoadRecord();
      m_dlbTimers[i].resetIdleRecord();
    }
  }

  /// \brief Return the load record of a DLB timer
  inline MFloat returnLoadRecord(const MInt dlbTimerId, const MInt mode = 0) {
    if(m_ignoreDlbTimers) {
      return -1.0;
    }
    ENSURE_VALID_TIMERID(dlbTimerId, name);
    ENSURE_VALID_TIMERMODE(mode, name);
    return m_dlbTimers[dlbTimerId].returnLoadRecord(mode);
  }

  /// \brief Return the idle record of a DLB timer
  inline MFloat returnIdleRecord(const MInt dlbTimerId, const MInt mode = 0) {
    if(m_ignoreDlbTimers) {
      return -1.0;
    }
    ENSURE_VALID_TIMERID(dlbTimerId, name);
    ENSURE_VALID_TIMERMODE(mode, name);
    return m_dlbTimers[dlbTimerId].returnIdleRecord(mode);
  }

  /// \brief Return the number of DLB timers
  MInt noDlbTimers() const { return m_dlbTimers.size(); }

  /// \brief Return the number of (sub-)timers for each DLB timer
  MInt noSubTimers() const { return maia::dlbTimer::DlbTimer::noTimers(); }

 private:
  /// Storage of DLB timers for all solvers/couplers/...
  std::vector<maia::dlbTimer::DlbTimer> m_dlbTimers{};
  /// Id of the currently running DLB load/idle timer
  MInt m_runningTimerId = -1;
  /// Current status of all timers; false: all timers disabled; true: at least one timer enabled
  MBool m_enabled = false;
  /// Global switch in createDlbTimers() to ignore all DLB timers, i.e. they cannot be enabled
  MBool m_ignoreDlbTimers = true;
};

#endif // DLBTIMER_H
