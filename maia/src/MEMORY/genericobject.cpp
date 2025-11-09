// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "genericobject.h"

#include <fstream>
#include <iomanip>
#include <typeinfo>
#include "COMM/mpioverride.h"
#include "IO/infoout.h"
#include "UTIL/debug.h"
#include "UTIL/timer.h"
#include "globalvariables.h"

using namespace std;
using namespace maia::alloc;

/** \brief Prints currently allocated memory
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
void printAllocatedMemory(const MLong oldAllocatedBytes, const MString& solverName, const MPI_Comm comm) {
  TRACE();

  m_log << fixed << "=================================================" << endl;
  m_log.precision(6);
  m_log << getMemorySize(g_allocatedBytes - oldAllocatedBytes) << " allocated by " << solverName << "." << endl;
  m_log << "Total memory: " << getMemorySize(g_allocatedBytes) << "." << endl;
  m_log << "=================================================" << endl;

  MLong allocBytes = g_allocatedBytes;
  MLong oldAllocBytes = oldAllocatedBytes;
  MLong allocMax = allocBytes - oldAllocatedBytes;
  if(globalNoDomains() > 1) {
    MPI_Allreduce(MPI_IN_PLACE, &allocBytes, 1, MPI_LONG, MPI_SUM, comm, AT_, "MPI_IN_PLACE", "allocBytes");
    MPI_Allreduce(MPI_IN_PLACE, &oldAllocBytes, 1, MPI_LONG, MPI_SUM, comm, AT_, "MPI_IN_PLACE", "oldAllocBytes");
    MPI_Allreduce(MPI_IN_PLACE, &allocMax, 1, MPI_LONG, MPI_MAX, comm, AT_, "MPI_IN_PLACE", "oldAllocBytes");
  }
  if(globalDomainId() == 0) {
    cerr.precision(6);
    cerr << "=== " << getMemorySize(allocBytes - oldAllocBytes) << " globally allocated by " << solverName
         << ". Total global memory: " << getMemorySize(allocBytes) << ". ===" << endl;
    cerr << "=== " << getMemorySize(allocMax) << " maximum globally allocated by " << solverName << ". ===" << endl;
  }
}

MString getMemorySize(MLong noBytes) {
  stringstream size;
  size.str("");
  char buffer[32];
  if((MFloat)noBytes / (1024 * 1024) < 1.0) {
    MFloat tmp = (MFloat)noBytes / (1024);
    sprintf(buffer, "%6.2f", tmp);
    size << buffer << " KB";
  } else if(((MFloat)noBytes / (1024 * 1024 * 1024)) > 1.0) {
    MFloat tmp = (MFloat)noBytes / (1024 * 1024 * 1024);
    sprintf(buffer, "%6.2f", tmp);
    size << buffer << " GB";
  } else {
    MFloat tmp = (MFloat)noBytes / (1024 * 1024);
    sprintf(buffer, "%6.2f", tmp);
    size << buffer << " MB";
  }
  return size.str();
}
