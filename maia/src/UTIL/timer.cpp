// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "timer.h"
#include "typetraits.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
//--------------------------------------------------------------------------------

MTimers& timers() {
  static MTimers timers;
  return timers;
}

vector<FunctionTiming> Profile::s_functionTimings;


MBool operator<(const FunctionTiming& a, const FunctionTiming& b) {
  if(approx(a.getDeltaCpuTime(), b.getDeltaCpuTime(), MFloatEps)) return (a.getInitCpuTime() < b.getInitCpuTime());
  return (a.getDeltaCpuTime() < b.getDeltaCpuTime());
}

Profile::~Profile() {
#ifdef _OPENMP
  if(omp_get_max_threads() > 1) {
    m_log << "Skipping profile output for OpenMP. FIXME: Tracer not thread-safe!" << std::endl;
    return;
  }
#endif

  const MFloat exitCpuTime = cpuTime();
  const MFloat exitWallTime = wallTime();
  const MFloat thresholdPercentage = 0.5;
  stringstream sstream;
  sstream << "    CPU      WALL   FUNCTION                    >> profile: '" << m_name << "' <<";
  const string header = sstream.str();
  const MString dashes1(header.size(), '_');
  const MString dashes2(header.size(), '-');
  m_log << dashes1 << endl << header << endl << dashes2 << endl << endl;

  MInt counter = 0;
  MInt supCounter = 0;
  // TODO labels:TIMERS add mode to compute statistics over all ranks
  if(s_functionTimings.size() > 0) {
    sort(s_functionTimings.begin(), s_functionTimings.end());
    reverse(s_functionTimings.begin(), s_functionTimings.end());
    for(vector<FunctionTiming>::size_type i = 0; i < s_functionTimings.size(); i++) {
      if(s_functionTimings[i].getInitCpuTime() < m_initCpuTime) continue;
      const MFloat relCpuTime =
          100.0 * s_functionTimings[i].getDeltaCpuTime() / max(1e-15, (exitCpuTime - m_initCpuTime));
      const MFloat relWallTime =
          100.0 * s_functionTimings[i].getDeltaWallTime() / max(1e-15, (exitWallTime - m_initWallTime));
      if(relCpuTime < thresholdPercentage) {
        supCounter++;
        continue;
      }
      char buffer[7];
      sprintf(buffer, "%6.2f", relCpuTime);
      char buffer2[7];
      sprintf(buffer2, "%6.2f", relWallTime);
      // TODO labels:TIMERS abbreviate function names/remove templates etc
      m_log << buffer << "%   " << buffer2 << "%   " << s_functionTimings[i].getName() << endl;
      counter++;
    }
    if(supCounter > 0) {
      m_log << "  .....     .....   (" << supCounter << " shorter timings with CPU<" << thresholdPercentage
            << "% were suppressed)" << endl;
    }
  }
  if(counter == 0) {
    m_log << "No timings recorded for timer '" << m_name << "'." << endl;
  }
  m_log << dashes2 << endl;
  m_log << "Total cpu time:  " << printTime(exitCpuTime - m_initCpuTime) << endl;
  m_log << "Total wall time: " << printTime(exitWallTime - m_initWallTime) << endl;
  m_log << dashes1 << endl;
}

MInt Profile::getTimingId(const MString name) {
  MInt tId = -1;
  if(s_functionTimings.size() > 0) {
    for(vector<FunctionTiming>::size_type i = 0; i < s_functionTimings.size(); i++) {
      if(s_functionTimings[i].getName() == name) {
        tId = i;
        break;
      }
    }
  }
  if(tId < 0) {
    tId = static_cast<MInt>(s_functionTimings.size());
    s_functionTimings.push_back(FunctionTiming(name));
  }
  ASSERT(tId > -1, "Non-existing timer");
  return tId;
}

MString Profile::printTime(MFloat secs) {
  stringstream time;
  time.str("");
  MFloat rem = secs;
  if(rem > 86400.0) {
    const MFloat div = floor(rem / 86400.0);
    time << ((MInt)div) << " days, ";
    rem -= div * 86400.0;
  }
  if(rem > 3600.0) {
    const MFloat div = floor(rem / 3600.0);
    time << ((MInt)div) << " hours, ";
    rem -= div * 3600.0;
  }
  if(rem > 60.0) {
    const MFloat div = floor(rem / 60.0);
    time << ((MInt)div) << " mins, ";
    rem -= div * 60.0;
  }
  time << rem << " secs";
  const MString ret = time.str();
  return ret;
}


FunctionTiming::FunctionTiming(string name)
  : m_initCpuTime(cpuTime()),
    m_deltaCpuTime(0),
    m_tmpCpuTime(0),
    m_initWallTime(wallTime()),
    m_deltaWallTime(0.0),
    m_tmpWallTime(-1.0),
    m_name(name) {}

FunctionTiming::~FunctionTiming() { m_name = "<deleted>"; }

FunctionTiming& FunctionTiming::operator=(const FunctionTiming& t) {
  m_initCpuTime = t.m_initCpuTime;
  m_deltaCpuTime = t.m_deltaCpuTime;
  m_tmpCpuTime = t.m_tmpCpuTime;
  m_initWallTime = t.m_initWallTime;
  m_deltaWallTime = t.m_deltaWallTime;
  m_tmpWallTime = t.m_tmpWallTime;
  m_name = t.m_name;
  return *this;
}

void FunctionTiming::in() {
  m_tmpCpuTime = cpuTime();
  m_tmpWallTime = wallTime();
}

void FunctionTiming::out() {
  if(m_tmpCpuTime > 0) m_deltaCpuTime += (cpuTime() - m_tmpCpuTime);
  if(m_tmpWallTime > 0.0) m_deltaWallTime += (wallTime() - m_tmpWallTime);
  m_tmpCpuTime = -1.0;
  m_tmpWallTime = -1.0;
}


/// \brief Output the min/max/average duration of a code section over the ranks in a communicator
/// Note: only use this function for timing initialization steps or similar due to the blocking MPI_Allreduce operations
/// involved. Instead use logDurations() to assemble a vector of durations first and evaluate at once.
///
/// \param[in] timeStart start time obtained with wallTime()/MPI_Wtime() to compute difference to
void logDuration_(const MFloat timeStart, const MString module, const MString comment, const MPI_Comm comm,
                  const MInt domainId, const MInt noDomains) {
  const MFloat duration = wallTime() - timeStart;

  std::vector<std::pair<MFloat, MString>> durations{};
  durations.push_back(std::make_pair(duration, comment));
  logDurations(durations, module, comm, domainId, noDomains);
}


/// \brief Output the min/max/average durations of provided timed code sections over the ranks in a communicator
void logDurations(std::vector<std::pair<MFloat, MString>>& durations, const MString module, const MPI_Comm comm,
                  const MInt domainId, const MInt noDomains) {
  const MInt noDurations = durations.size();

  std::vector<MFloat> maxDurations(noDurations);
  for(MInt i = 0; i < noDurations; i++) {
    maxDurations[i] = durations[i].first;
  }
  // Copy durations vector
  std::vector<MFloat> minDurations = maxDurations;
  std::vector<MFloat> sumDurations = maxDurations;

  // Compute max, min and sum of durations over all involved ranks
  MPI_Allreduce(MPI_IN_PLACE, &maxDurations[0], noDurations, maia::type_traits<MFloat>::mpiType(), MPI_MAX, comm, AT_,
                "MPI_IN_PLACE", "maxDurations");
  MPI_Allreduce(MPI_IN_PLACE, &minDurations[0], noDurations, maia::type_traits<MFloat>::mpiType(), MPI_MIN, comm, AT_,
                "MPI_IN_PLACE", "minDurations");
  MPI_Allreduce(MPI_IN_PLACE, &sumDurations[0], noDurations, maia::type_traits<MFloat>::mpiType(), MPI_SUM, comm, AT_,
                "MPI_IN_PLACE", "sumDurations");

  const MInt maxLineLength = 256;
  MChar b[maxLineLength];
  for(MInt i = 0; i < noDurations; i++) {
    const MString comment = durations[i].second;
    snprintf(b, maxLineLength, "=== MAIA %s DURATION: %-35s | min: %.4e s | avg: %.4e s | max: %.4e s |",
             module.c_str(), comment.c_str(), minDurations[i], sumDurations[i] / (MFloat)noDomains, maxDurations[i]);
    if(domainId == 0) {
      std::cerr << b << std::endl;
    }
    m_log << b << std::endl;
  }
}
