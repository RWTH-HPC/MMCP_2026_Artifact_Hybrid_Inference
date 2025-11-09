// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "functions.h"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include "backtrace.h"
#include "environment.h"
#include "globals.h"

#ifndef PVPLUGIN
extern Environment* mEnvironment;
#endif

/** Terminates MAIA properly
 * @author Pascal Meysonnat
 *
 */

void mTerm(const MInt errorCode, const MString& location, const MString& message) {
  if(errorCode != 0) {
    std::stringstream s;
    s << "\n";
    s << "Rank " << globalDomainId() << " threw exit code " << errorCode << "\n";
    s << "Error in " << location << ": " << message << "\n";
    s << "\n"
      << "Program is aborting!!\n";
    std::cerr << s.str() << std::flush;

    // Print backtrace (if enabled)
    BACKTRACE();

    // Close the log file to make sure that no MPI error occurs from the
    // unclosed file, and that a proper XML footer is written
#ifndef PVPLUGIN
    m_log.close(true);
    maia_res.close(true);
#endif
    MPI_Abort(globalMaiaCommWorld(), errorCode, AT_);
  } else {
    mDealloc();
#ifndef PVPLUGIN
    if(mEnvironment) {
      delete mEnvironment;
      mEnvironment = nullptr;
    }

    // Close the log file to make sure that no MPI error occurs from the
    // unclosed file, and that a proper XML footer is written
    m_log.close();
    maia_res.close();
#endif
    // Call MPI_Finalize to ensure proper MPI shutdown
    MPI_Finalize();

    // Exit the program
    exit(0);
  }
  exit(errorCode);
}


/// \brief Returns true if the file \p fileName exists, false otherwise.
MBool fileExists(const MString& fileName) {
  struct stat buffer;
  return (stat(fileName.c_str(), &buffer) == 0);
}

/// \brief Copies file \p fromName to file \p toName
///
/// \returns error code
/// \retval 0 if copying succeeded
/// \retval -1 if copying failed
MInt copyFile(const MString& fromName, const MString& toName) {
  if(!fileExists(fromName)) {
    std::cerr << "Could not copy file " << fromName << ".\n";
    return -1;
  }
  if(fileExists(toName)) {
    remove(toName.c_str());
  }
  std::ifstream src(fromName.c_str());
  std::ofstream dst(toName.c_str());
  if(src.good() && dst.good()) {
    dst << src.rdbuf();
    dst.close();
    dst.clear();
    src.close();
    src.clear();
  } else {
    std::cerr << "Could not copy file " << fromName << " (2).\n";
    return -1;
  }
  return 0;
}


/// \brief Checks if the given grid extents and cell sizes match when creating a multisolver grid and
///        corresponding single solver grids to ensure that the same min-level cells and Hilbert
///        order is used
///
/// Note: from createGridMap()
void checkMultiSolverGridExtents(const MInt nDim, const MFloat* centerOfGravity, const MFloat lengthLevel0,
                                 const MInt minLevel, const MFloat* targetGridCenterOfGravity,
                                 const MFloat targetGridLengthLevel0, const MInt targetGridMinLevel) {
  // Define epsilon for floating point comparisons (argh!)... the original
  // definition is taken from the grid constructor but it must be clear to anyone
  // that this is a less-than-optimal solution
  const MFloat eps = 1.0 / FPOW2(30) * lengthLevel0;

  m_log << std::setprecision(15) << "Check multisolver grid extends: length0=" << lengthLevel0
        << "; minLevel=" << minLevel << "; targetLenght0=" << targetGridLengthLevel0
        << "; targetMinLevel=" << targetGridMinLevel << "; eps=" << eps << std::endl;

  // Check extents
  const std::array<MString, 3> dirs = {{"x", "y", "z"}};
  for(MInt dir = 0; dir < nDim; dir++) {
    const std::array<MFloat, 2> gridExtent = {
        {centerOfGravity[dir] - 0.5 * lengthLevel0, centerOfGravity[dir] + 0.5 * lengthLevel0}};
    const std::array<MFloat, 2> globalExtent = {{targetGridCenterOfGravity[dir] - 0.5 * targetGridLengthLevel0,
                                                 targetGridCenterOfGravity[dir] + 0.5 * targetGridLengthLevel0}};

    if(gridExtent[0] + eps < globalExtent[0]) {
      TERMM(1,
            "Grid extents exceed multisolver bouding box in negative " + dirs[dir]
                + "-direction: grid=" + std::to_string(gridExtent[0]) + "; box=" + std::to_string(globalExtent[0]));
    }
    if(gridExtent[1] - eps > globalExtent[1]) {
      TERMM(1,
            "Grid extents exceed multisolver bouding box in positive " + dirs[dir]
                + "-direction: grid=" + std::to_string(gridExtent[1]) + "; box=" + std::to_string(globalExtent[1]));
    }
    m_log << std::setprecision(15) << "Grid extents in the " << dirs[dir] << "-direction: [" << gridExtent[0] << ", "
          << gridExtent[1] << "]; global: [" << globalExtent[0] << ", " << globalExtent[1] << "]" << std::endl;
  }

  // Check cell length at min level
  const MFloat gridMinLevelLength = lengthLevel0 * FFPOW2(minLevel);
  const MFloat globalMinLevelLength = targetGridLengthLevel0 * FFPOW2(targetGridMinLevel);
  if(fabs(gridMinLevelLength - globalMinLevelLength) > eps) {
    TERMM(1,
          "Length of min level cells do not match between grid and given multisolver bounding "
          "box: grid: "
              + std::to_string(gridMinLevelLength) + "; box: " + std::to_string(globalMinLevelLength));
  } else {
    m_log << std::setprecision(15)
          << "Length of min level cells match between grid and given multisolver bounding "
             "box: grid: "
          << gridMinLevelLength << "; box: " << globalMinLevelLength << std::endl;
  }

  // Check if grid centers are displaced by an integer multiple of the min
  // level length
  for(MInt dir = 0; dir < nDim; dir++) {
    const MFloat displacement = fabs(centerOfGravity[dir] - targetGridCenterOfGravity[dir]);
    const MFloat quotient = displacement / globalMinLevelLength;
    if(!isApproxInt(quotient, eps)) {
      TERMM(1, "The grid centers are displaced in the " + dirs[dir]
                   + "-direction by a non-integer multiple of the length of a "
                     "partition cell: "
                   + std::to_string(quotient));
    } else {
      m_log << std::setprecision(15) << "The grid centers are displaced in the " << dirs[dir]
            << "-direction by a multiple of the length of a "
               "partition cell: "
            << quotient << " (displacement = " << displacement << ")" << std::endl;
    }
  }
}


/// \brief Loads point coordinates from an input file.
MInt loadPointCoordinatesFromFile(const MString inputFileName, const MInt nDim, std::vector<MFloat>& coordinates) {
  TRACE();

  MString line;
  MFloat curFloat;
  MInt noPoints = 0;
  std::istringstream iss;
  std::ifstream csvFile(inputFileName);
  coordinates.clear();

  // Read all lines and get coordinates
  while(getline(csvFile, line)) {
    iss.str(line);
    iss.clear();

    for(MInt i = 0; i < nDim; i++) {
      iss >> curFloat;
      if(iss.fail()) {
        std::ostringstream err;
        // start line count at one
        err << "Error at line " << noPoints + 1 << ": " << line << "\n"
            << "Either wrong dimension (nDim = " << nDim << ") or otherwise wrong format."
            << "Format should be nDim floats seperated by spaces per line.";
        TERMM(1, err.str());
      }
      coordinates.push_back(curFloat);
    }
    noPoints++;
  }

  return noPoints;
}


#ifndef PVPLUGIN
/// \brief Write memory statistics
void writeMemoryStatistics(const MPI_Comm comm, const MInt noDomains, const MInt domainId, const MString at,
                           const MString comment) {
  // Static members to compute difference to last memory evaluation
  static MLong virtMemLast = 0;
  static MLong physMemLast = 0;
  // OBTAIN MEMORY USAGE
  MLong virtMem = 0;
  MLong physMem = 0;
  MLong stackMem = 0;
  MLong physMemFree = 0;
  MLong memAvailable = 0;

  {
    // Create error flags
    MInt fileNotFound = 0;
    MInt memoryNotFound = 0;

    // Open status file of the process
    std::ifstream fin;
    fin.open("/proc/self/status");
    // Modify flag if file is not found
    if(!fin) {
      fileNotFound = 1;
    }

    MString line;
    MString name;
    MLong buffer;
    std::istringstream iss;
    std::array<MInt, 5> foundInfo;
    foundInfo.fill(false);

    // Read all lines and get memory usage/allocation size for this process
    while(getline(fin, line)) {
      buffer = 0;
      iss.str(line);
      iss.clear();
      getline(iss, name, ':');
      if(name == "VmRSS") { // Physical memory currently used by current process
        iss >> buffer;
        physMem = buffer;
        foundInfo[0] = true;
      } else if(name == "VmData") { // Allocated memory size on heap
        iss >> buffer;
        virtMem += buffer;
        foundInfo[1] = true;
      } else if(name == "VmStk") { // Allocated memory size on stack
        iss >> buffer;
        virtMem += buffer;
        stackMem = buffer;
        foundInfo[2] = true;
      }
    }
    fin.close();

    // Open meminfo file on current node to determine total free/available memory per node
    fin.open("/proc/meminfo");
    if(!fin) {
      fileNotFound = 1;
    }

    while(getline(fin, line)) {
      buffer = 0;
      iss.str(line);
      iss.clear();
      getline(iss, name, ':');
      if(name == "MemFree") {
        iss >> buffer;
        physMemFree = buffer;
        foundInfo[3] = true;
      } else if(name == "MemAvailable") {
        iss >> buffer;
        memAvailable = buffer;
        foundInfo[4] = true;
      }
    }
    fin.close();

    // Set flag if any memory usage couldn't be gathered from file
    memoryNotFound = std::any_of(foundInfo.begin(), foundInfo.end(), [](MBool i) { return !i; });

    // Allreduce flags, only proceed if there was no error
    MPI_Allreduce(MPI_IN_PLACE, &fileNotFound, 1, MPI_INT, MPI_MAX, comm, AT_, "MPI_IN_PLACE", "&fileNotFound");
    MPI_Allreduce(MPI_IN_PLACE, &memoryNotFound, 1, MPI_INT, MPI_MAX, comm, AT_, "MPI_IN_PLACE", "&memoryNotFound");
    MPI_Allreduce(MPI_IN_PLACE, &foundInfo[0], 5, MPI_INT, MPI_MAX, comm, AT_, "MPI_IN_PLACE", "&foundInfo");

    if(fileNotFound || memoryNotFound) {
      // Throw error and return from function if file or memory couldn't be accessed
      if(domainId == 0) {
        std::stringstream ss;
        ss << "Error in writeMemoryStatistics: Could not determine memory statistics! " << fileNotFound << " "
           << memoryNotFound << " (";
        for(MInt i = 0; i < 5; i++) {
          ss << " " << foundInfo[i];
        }
        ss << ")" << std::endl;
        std::cerr << ss.str();
        m_log << ss.str();
      }
      return;
    }
  }

  MLongScratchSpace physMemPerProcess(noDomains, AT_, "physMemPerProcess");
  MLongScratchSpace virtMemPerProcess(noDomains, AT_, "virtMemPerProcess");

  // Gather memory from each process
  MPI_Gather(&physMem, 1, MPI_LONG, physMemPerProcess.getPointer(), 1, MPI_LONG, 0, comm, AT_, "physMem",
             "physMemPerProcess.getPointer()");
  MPI_Gather(&virtMem, 1, MPI_LONG, virtMemPerProcess.getPointer(), 1, MPI_LONG, 0, comm, AT_, "virtMem",
             "virtMemPerProcess.getPointer()");

  MLongScratchSpace physMemDiffPerProcess(noDomains, AT_, "physMemDiffPerProcess");
  MLongScratchSpace virtMemDiffPerProcess(noDomains, AT_, "virtMemDiffPerProcess");

  MLong physMemDiff = physMem - physMemLast;
  MLong virtMemDiff = virtMem - virtMemLast;
  // Gather difference in allocated memory from each process
  MPI_Gather(&physMemDiff, 1, MPI_LONG, physMemDiffPerProcess.getPointer(), 1, MPI_LONG, 0, comm, AT_, "physMem",
             "physMemDiffPerProcess.getPointer()");
  MPI_Gather(&virtMemDiff, 1, MPI_LONG, virtMemDiffPerProcess.getPointer(), 1, MPI_LONG, 0, comm, AT_, "virtMem",
             "virtMemDiffPerProcess.getPointer()");

  // Store current values for next evaluation
  physMemLast = physMem;
  virtMemLast = virtMem;

  // Gather stack memory usage
  MLongScratchSpace stackMemPerProcess(noDomains, AT_, "stackMemPerProcess");
  MPI_Gather(&stackMem, 1, MPI_LONG, stackMemPerProcess.getPointer(), 1, MPI_LONG, 0, comm, AT_, "stackMem",
             "stackMemPerProcess.getPointer()");

  // Compute minimum free/available memory
  MLong minPhysMemFree = physMemFree;
  MPI_Allreduce(MPI_IN_PLACE, &minPhysMemFree, 1, MPI_LONG, MPI_MIN, comm, AT_, "MPI_IN_PLACE", "&minPhysMemFree");
  MLong minMemAvailable = memAvailable;
  MPI_Allreduce(MPI_IN_PLACE, &minMemAvailable, 1, MPI_LONG, MPI_MIN, comm, AT_, "MPI_IN_PLACE", "&minMemAvailable");

  // Determine global memory statistics
  if(domainId == 0) {
    MLong totalPhysMem = 0;
    MLong totalVirtMem = 0;
    MLong totalPhysMemDiffSum = 0;
    MLong totalVirtMemDiffSum = 0;
    MLong totalPhysMemDiffMax = 0;
    MLong totalVirtMemDiffMax = 0;
    MLong maxStackMem = 0;

    for(MInt i = 0; i < noDomains; i++) {
      // Total allocation size
      totalPhysMem += physMemPerProcess[i];
      totalVirtMem += virtMemPerProcess[i];
      // Total difference in allocation size compared to the last evaluation
      totalPhysMemDiffSum += physMemDiffPerProcess[i];
      totalVirtMemDiffSum += virtMemDiffPerProcess[i];
      // Maximum difference in allocation size over all processes
      totalPhysMemDiffMax = std::max(physMemDiffPerProcess[i], totalPhysMemDiffMax);
      totalVirtMemDiffMax = std::max(virtMemDiffPerProcess[i], totalVirtMemDiffMax);
      // Maximum stack memory usage
      maxStackMem = std::max(stackMemPerProcess[i], maxStackMem);
    }

    // Get maximum size of the process stack (this assumes it is the same on all ranks)
    rlimit rlim;
    getrlimit(RLIMIT_STACK, &rlim);
    // Define rlimit as an unsigned long (same as it's defined in the kernel)
    const MUlong rlim_stack = rlim.rlim_cur;
    // Note: issue warning if stack usage is quite high (compared to stack size limit)
    if(maxStackMem > 0.5 * rlim_stack) {
      std::stringstream warning;
      warning << std::endl
              << "WARNING: maximum stack memory usage >50% of its limit, use 'ulimit -s unlimited' to remove this "
                 "memory restriction and avoid segmentation faults if the stack memory usage exceeds its limit."
              << std::endl;
      warning << "WARNING: stack memory usage " << (MFloat)maxStackMem << " KB; stack limit "
              << (MFloat)rlim_stack / 1024 << " KB" << std::endl
              << std::endl;
      m_log << warning.str();
      std::cerr << warning.str();
    }

    // Write memory statistics
#ifndef NDEBUG
    // Memory per process
    for(MInt i = 0; i < noDomains; i++) {
      m_log << " Process " << i << " - Current memory usage: physical = " << (MFloat)physMemPerProcess[i] / 1024
            << " MB; allocation = " << (MFloat)virtMemPerProcess[i] / 1024 << " MB" << std::endl;
    }
#endif

    std::stringstream ss;
    ss << std::endl;
    ss << "******************************* MEMORY STATISTICS *******************************" << std::endl;
    ss << "***** Comment: " << comment << " - #ranks: " << noDomains << std::endl;
    ss << "***** Location: " << at << std::endl;
    ss << "***** " << std::endl;
    // Average memory
    ss << "***** Average memory usage: physical = " << (MFloat)totalPhysMem / (noDomains * 1024)
       << " MB; allocation = " << (MFloat)totalVirtMem / (noDomains * 1024) << " MB\n";

    // Min/Max memory
    ss << "***** Minimun memory usage: physical = "
       << (MFloat)*std::min_element(physMemPerProcess.begin(), physMemPerProcess.end()) / 1024
       << " MB; allocation = " << (MFloat)*std::min_element(virtMemPerProcess.begin(), virtMemPerProcess.end()) / 1024
       << " MB\n";
    ss << "***** Maximum memory usage: physical = "
       << (MFloat)*std::max_element(physMemPerProcess.begin(), physMemPerProcess.end()) / 1024
       << " MB; allocation = " << (MFloat)*std::max_element(virtMemPerProcess.begin(), virtMemPerProcess.end()) / 1024
       << " MB\n";
    ss << "***** Maximum diff in memory usage: physical = " << (MFloat)totalPhysMemDiffMax / 1024
       << " MB; allocation = " << (MFloat)totalVirtMemDiffMax / 1024 << " MB\n";

    // Total memory
    ss << "***** Total physical memory usage (RAM): " << (MFloat)totalPhysMem / (1024 * 1024) << " GB\n";
    ss << "***** Diff total physical memory usage (RAM): " << (MFloat)totalPhysMemDiffSum / (1024 * 1024) << " GB\n";
    ss << "***** Total allocation size (Virtual Memory): " << (MFloat)totalVirtMem / (1024 * 1024) << " GB"
       << std::endl;
    ss << "***** Diff total allocation size (Virtual Memory): " << (MFloat)totalVirtMemDiffSum / (1024 * 1024) << " GB"
       << std::endl;
    ss << "***** " << std::endl;
    ss << "***** Maximum stack memory: " << (MFloat)maxStackMem << " KB; stack limit " << (MFloat)rlim_stack / 1024
       << " KB" << std::endl;
    ss << "***** " << std::endl;
    ss << "***** Minimum available memory per node (meminfo): " << (MFloat)minMemAvailable / (1024 * 1024) << " GB"
       << std::endl;
    ss << "***** Minimum free memory per node (RAM): " << (MFloat)minPhysMemFree / (1024 * 1024) << " GB" << std::endl;
    ss << "******************************* MEMORY STATISTICS *******************************" << std::endl << std::endl;

    std::cout << ss.str() << std::endl;
    m_log << ss.str() << std::endl;
  }
}
#endif
