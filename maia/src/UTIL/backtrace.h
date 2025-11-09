// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef BACKTRACE_H_
#define BACKTRACE_H_

#include <iostream>
#include <sstream>
#include "INCLUDE/maiatypes.h"
#include "compiler_config.h"

#ifdef MAIA_GCC_COMPILER
// Needed for stack trace
#include <cxxabi.h>
#include <execinfo.h>
#endif

// General backtrace macro
#if defined(ENABLE_BACKTRACE)
#include "llvm/Support/Signals.h"
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include "llvm/Support/raw_ostream.h"

#define BACKTRACE()                                                                                                    \
  do {                                                                                                                 \
    maia::debug::backtrace();                                                                                          \
  } while(false)
#else
#define BACKTRACE()                                                                                                    \
  do {                                                                                                                 \
  } while(false)
#endif


namespace maia {
namespace debug {

/**
 * \brief Prints a backtrace of the function call path if possible.
 *
 * \author Fabian Klemp <f.klemp@aia.rwth-aachen.de>
 * \date 2015-09-16
 *
 * \details Uses the LLVMSupport library to print a stacktrace
 *          Note: use BACKTRACE(...) instead of calling this method directly
 *
 */

#if defined(ENABLE_BACKTRACE)
inline void backtrace() {
  llvm::errs() << "Backtrace (line numbers may be too large by 1-3 lines):\n";
  llvm::sys::PrintStackTrace(llvm::errs());
}
#endif

/**
 * \brief Prints a backtrace of the function call path if possible.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-10-09
 *
 * \param[in] noFramesOmitted Set number of omitted frames (e.g. to hide
 *                            backtrace() itself).
 *
 * \details Works only with GCC. Needs command line tool 'addr2line' to be
 *          callable from within MAIA.
 */
#ifdef MAIA_GCC_COMPILER
inline void backtrace_gcc(const MString& fileName, const MInt noFramesOmitted = 1) {
  // Get trace
  const int maxNoTraces = 128;
  void* trace[maxNoTraces];
  int noTraces = ::backtrace(trace, maxNoTraces);

  // Get messages (i.e. function names)
  char** messages = backtrace_symbols(trace, noTraces);

  // Create stream and char buffer for single-flush printing
  std::stringstream s;
  const MInt maxLineLength = 1024;
  MChar line[maxLineLength];

  // Print out messages one by one
  const MInt maxCmdSize = 1024;
  MInt frameId = 0;
  const MInt noFrames = noTraces - noFramesOmitted;
  const MInt width = noFrames ? floor(log10(noFrames)) + 1 : 1;
  s << "\nBacktrace (line numbers are usually too large, typically by 1-3 "
       "lines):\n";
  for(MInt i = noFramesOmitted; i < noTraces; i++) {
    // Get demangled function name
    MString mangled(messages[i]);
    mangled = mangled.substr(mangled.find("(") + 1);
    mangled = mangled.substr(0, mangled.find("+"));
    int status;
    char* demangledChar = abi::__cxa_demangle(mangled.c_str(), 0, 0, &status);

    // Only use demangled name if demangling worked
    const MString demangled = (status == 0) ? demangledChar : messages[i];
    free(demangledChar);

    // Get filename/line number using addr2line
    char cmd[maxCmdSize];
    std::sprintf(cmd, "addr2line %p -s -e %s", trace[i], fileName.c_str());
    MString fileline;
    FILE* pipe = popen(cmd, "r");
    if(!pipe) {
      fileline = "n/a";
    } else {
      char buffer[128];
      while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != nullptr) {
          fileline += buffer;
        }
      }
    }
    pclose(pipe);

    // Parse output of addr2line to skip newline and possible discriminator
    // information
    fileline = fileline.substr(0, fileline.length() - 1);
    fileline = fileline.substr(0, fileline.find(" "));

    // Print function names + filename:lineno
    std::sprintf(line, "#%-*d 0x%016lx in %s at %s\n", width, frameId++, reinterpret_cast<uintptr_t>(trace[i]),
                 demangled.c_str(), fileline.c_str());
    line[maxLineLength - 1] = '\0';
    s << line;
  }
  std::cout << s.str() << std::endl;
}
#else
inline void backtrace_gcc(const MString& NotUsed(fileName), MInt NotUsed(noFramesOmitted)) {
  std::cout << "Not using GCC - backtrace disabled." << std::endl;
}
#endif

} // namespace debug
} // namespace maia

#endif // ifndef BACKTRACE_H_
