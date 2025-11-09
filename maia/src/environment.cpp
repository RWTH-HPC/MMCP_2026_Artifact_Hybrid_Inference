// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "environment.h"

#include <array>
#if defined(MAIA_MS_COMPILER)
#include <Windows.h>
#include <Winsock2.h>
#include <direct.h>
#else
#include <dirent.h>
#ifndef _SX
#include <getopt.h>
#endif
#include <pwd.h>
#include <unistd.h>
#endif
#include <ctime> // Needed for generating time strings
#include <fcntl.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "application.h"
#include "globals.h"

#ifdef _SX
#include <sys/socket.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

//! Reads the name of the property-file and creates a new Application

MInt Environment::m_argc;
MChar** Environment::m_argv;

Environment::Environment(MInt argc, char** argv) {
  TRACE();

  NEW_TIMER_GROUP(tg_newenv, "New environment");
  NEW_TIMER(t_completenewenv, "complete new environment", tg_newenv);
  RECORD_TIMER_START(t_completenewenv);

  // Copy command line parameters to member variables and process them
  NEW_SUB_TIMER(t_parseCL, "parse command line", t_completenewenv);
  RECORD_TIMER_START(t_parseCL);
  Environment::m_argc = argc;
  Environment::m_argv = argv;
  parseCommandline();
  RECORD_TIMER_STOP(t_parseCL);

  // After parsing optional command line arguments, print startup information
  if(globalDomainId() == 0) {
    printStartupInformation();
  }

  // In order to read the property file via ParallelIo a (small) scratch is necessary
  // Since the properties which define the actual scratch are not known yet
  // a size of 1000 should be sufficient to perform all property related IO
  NEW_SUB_TIMER(t_readProp, "read property file", t_completenewenv);
  RECORD_TIMER_START(t_readProp);
  MInt tmpScratchNoCells = 750000;
  mScratch = new Scratch(1.0, tmpScratchNoCells);

  // Read in property file so that all properties are available from this point on
  DEBUG("Environment:: Context::readPropertyFile", MAIA_DEBUG_IO);
  if(m_propertyFileInput.find(".toml") == MString::npos) {
    Context::readPropertyFile(NETCDF, m_propertyFileInput);
  } else {
    Context::readPropertyFile(TOML, m_propertyFileInput);
  }
  // Delete the temporary scratch
  delete mScratch;

  RECORD_TIMER_STOP(t_readProp);

  // Read in properties on maximum number of cells to know how much scratch space needs to be allocated
  /*! \page propertyPage1
      \section maxNoCells
      <code>MInt maxNoCells</code> \n
      default = <code>nullptr</code> \n \n
      How many cells are maximal allowed, memory issue! \n
      Keywords: <i>GENERAL</i>
  */

#ifdef _OPENMP
  // use value from properties file, if not given, use value from environment, if not set, use omp=1
  MInt maxNoOMPThreads = 1;
  if(getenv("OMP_NUM_THREADS")) {
    maxNoOMPThreads = atoi(getenv("OMP_NUM_THREADS"));
  }
  maxNoOMPThreads = Context::getBasicProperty<MInt>("numOMPThreads", AT_, &maxNoOMPThreads);
  omp_set_num_threads(maxNoOMPThreads);
#endif

  MInt maxNoCells = Context::getBasicProperty<MInt>("maxNoCells", AT_);
  if(Context::propertyExists("noDomains")) {
    const MInt testNoDomains = Context::getBasicProperty<MInt>("noDomains", AT_, 0);
    if(globalNoDomains() < testNoDomains) {
      // Here, the number of maxNoCells is scaled. This is useful if a test
      // case is specified to run with a certain number of ranks
      // ('noDomains' property in run.toml) but needs to be run on a lower
      // number of mpiranks, p.e. by running maia on an accelerator.
      maxNoCells *= testNoDomains / (MFloat)globalNoDomains();
      cerr0 << "noDomain > number of used mpi ranks! Therefore, increasing maxNoCells: " << maxNoCells << std::endl;
    }
  }

  /*! \page propertyPage1
      \section scratchSize
      <code>MInt Scratch::m_number_of_elements</code> \n
      default = <code>1.0</code> \n \n
      Defines size of the scratch space in bytes: m_number_of_elements = scratchSize*maxNoCells*sizeOf(MFloat) bytes
     \n Keywords: <i>GENERAL</i>
  */

  NEW_SUB_TIMER(t_allocScratch, "allocate scratch", t_completenewenv);
  RECORD_TIMER_START(t_allocScratch);
  MFloat scratchSize = 2.0;
  scratchSize = Context::getBasicProperty<MFloat>("scratchSize", AT_, &scratchSize);

  // Allocate new scratch space
  mScratch = new Scratch(scratchSize, maxNoCells);
  m_log << Scratch::printSelf();
  RECORD_TIMER_STOP(t_allocScratch);

  // Create a new MAIA application, the class that controls the actual grid generator/flow solver
  DEBUG("Environment:: create Application", MAIA_DEBUG_LEVEL1);
  mApplication = new Application();
  RECORD_TIMER_STOP(t_completenewenv);
}

//! Destructor
Environment::~Environment() {
  m_log << Scratch::printSelf();
  delete mScratch;
  if(mApplication != nullptr) delete mApplication;
}

//! Runs the Environment, makes the Application run
MInt Environment::run() {
  Profile runProfile("Environment::run");

  Context::communicateProperties();
  Context::initializationProcessFinished();

  TRACE();

  MBool flowSolver = false;
  flowSolver = Context::getBasicProperty<MBool>("flowSolver", AT_, &flowSolver);

  // Return here if flow solver is not enabled.
  if(!flowSolver) {
    return 0;
  }

  MInt nDim = read_nDim();

  if(nDim == 2) {
    mApplication->run<2>();
  } else if(nDim == 3) {
    mApplication->run<3>();
  } else {
    mTerm(1, AT_, "Invalid value of nDim!");
  }

  return 0;
}

//! Ends the Environment
MInt Environment::end() {
  TRACE();

#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
  Context::writePropertiesHumanReadable();
#endif
  return 0;
}


//! Reads out the options that where used by the program call
void Environment::parseCommandline() {
  TRACE();

#if defined(MAIA_MS_COMPILER) || defined(_SX)
#pragma message("WARNING: Not implemented!")
#else
  opterr = 0;
  while(true) {
    int c = 0;
    DIR* directory = nullptr;
    string inputDirectory;
    string outputDirectory;
    int option_index = 0;
    static option long_options[] = {
        {"debug", 1, 0, 'd'},   {"help", 0, 0, 'h'},     {"input", 1, 0, 'i'},      {"output", 1, 0, 'o'},
        {"version", 0, 0, 'v'}, {"compiler", 0, 0, 'c'}, {"build-type", 0, 0, 'b'}, {0, 0, 0, 0}};

    c = getopt_long(m_argc, m_argv, "d:hi:o:p:vcb", long_options, &option_index);
    if(c == -1) {
      break;
    }

    switch(c) {
      case 0:
        cout << "option " << long_options[option_index].name << endl;
        if(optarg != nullptr) cout << "with arg " << optarg << endl;
        cout << endl;
        break;

      case 'd':
        if(atoi(optarg) < 2 * MDebug::m_maxLevel && atoi(optarg) >= MDebug::m_minLevel) {
          SET_DEBUG_LEVEL((MDebugLevel)atoi(optarg));
          if(globalDomainId() == 0) {
            cout << "SET DEBUG LEVEL TO: " << atoi(optarg) << endl;
          }
          m_log << "SET DEBUG LEVEL TO: " << atoi(optarg) << endl;
          DEBUG("Environment::parseCommandline: option -d " << atoi(optarg), MAIA_DEBUG_IO);
        } else {
          cerr << "Error: debug level out of range!" << endl;
          cerr << "Specified debug level is " << atoi(optarg) << ", but must be in the range of " << endl;
          cerr << 2 * MDebug::m_minLevel << " <= debug level <= " << 2 * MDebug::m_maxLevel << endl;
          cerr << "Debug level is unchanged" << endl;
          return;
        }
        break;

      case 'h':
        cout << "usage: " << m_argv[0] << " [-h] [-b] [-c] [--debug LEVEL] [--input DIRECTORY] "
             << "[--output DIRECTORY] [-v] [PROPERTY_FILE]" << endl
             << endl
             << "positional arguments:" << endl
             << "  PROPERTY_FILE            name of property file to use (default: '" << m_propertyFileInput << "')\n"
             << endl
             << "optional arguments:\n"
             << "  -b, --build-type         show build type that was used\n"
             << "  -c, --compiler           show used compiler\n"
             << "  -d, --debug LEVEL        use specified debug level\n"
             << "  -h, --help               this help screen\n"
             << "  -i, --input DIRECTORY    specified directory is used for all "
                "input files\n"
             << "  -o, --output DIRECTORY   specified directory is used for all "
                "output files\n"
             << "  -v, --version            show MAIA version\n"
             << endl;
        m_log.close();
        maia_res.close();
        MPI_Finalize();
        exit(0);
        break;

      case 'i':
        directory = opendir(optarg);
        if(directory != nullptr) {
          closedir(directory);
          inputDirectory.append(optarg);
          DEBUG("Environment::parseCommandline: option -i " << optarg, MAIA_DEBUG_IO);
        } else {
          stringstream errorMessage;
          errorMessage << " error opening inputDirectory \"" << optarg << "\"!" << endl;
          mTerm(1, AT_, errorMessage.str());
        }
        break;

      case 'o':
        directory = opendir(optarg);
        if(directory != nullptr) {
          closedir(directory);
          outputDirectory.append(optarg);
          DEBUG("Environment::parseCommandline: option -o " << optarg, MAIA_DEBUG_IO);
        } else {
          stringstream errorMessage;
          errorMessage << "Error opening outputDirectory \"" << optarg << "\"!";
          mTerm(1, AT_, errorMessage.str());
        }
        break;

      case 'v':
// Show version and quit
#ifndef MAIA_VERSION_STRING
#define MAIA_VERSION_STRING unknown
#endif
        cout << XSTRINGIFY(MAIA_VERSION_STRING) << endl;
        m_log.close();
        maia_res.close();
        MPI_Finalize();
        exit(0);
        break;

      case 'c':
// Show used compiler and quit
#ifndef MAIA_COMPILER_STRING
#define MAIA_COMPILER_STRING unknown
#endif
        cout << XSTRINGIFY(MAIA_COMPILER_STRING) << endl;
        m_log.close();
        maia_res.close();
        MPI_Finalize();
        exit(0);
        break;

      case 'b':
// Show used build type and quit
#ifndef MAIA_BUILD_TYPE_STRING
#define MAIA_BUILD_TYPE_STRING unknown
#endif
        cout << XSTRINGIFY(MAIA_BUILD_TYPE_STRING) << endl;
        m_log.close();
        maia_res.close();
        MPI_Finalize();
        exit(0);
        break;

      case '?':
        break;

      default:
        cerr << "Environment::parseCommandline: Error getopt returned "
                "unkown option code \""
             << c << "\"" << endl;
    }
  }

  // Obtain command line argument for property file
  if(optind < m_argc) {
    m_propertyFileInput = m_argv[optind];
  }
  // On domain 0, check if file is readable
  if(globalDomainId() == 0 && !ifstream(m_propertyFileInput)) {
    TERMM(1, "property file '" + m_propertyFileInput + "' not readable");
  }
#endif
}


/// \brief Print startup summary of MAIA.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-09-19
///
/// This method prints some information to stdout that may be useful for
/// debugging purposes later (and to identify which stdout logs belong to which
/// MAIA run).
void Environment::printStartupInformation() {
  // Define buffer size
  const MInt bufferSize = 1024;

  // Get current hostname
  array<char, bufferSize> host_{};
  gethostname(&host_[0], bufferSize - 1);
  host_[bufferSize - 1] = '\0';
  const MString host(&host_[0]);

  // Get current username
  MString user = "(unknown user)";
#if defined(MAIA_MS_COMPILER)
  constexpr MInt INFO_BUFFER_SIZE = 32767;
  TCHAR infoBuf[INFO_BUFFER_SIZE];
  DWORD bufCharCount = INFO_BUFFER_SIZE;
  if(GetUserName(infoBuf, &bufCharCount)) {
    user = infoBuf;
  }
#else
  passwd* p;
  p = getpwuid(getuid());
  if(p) {
    user = MString(p->pw_name);
  }
#endif

  // Get current directory
  array<char, bufferSize> dir_{};
#if defined(MAIA_MS_COMPILER)
  _getcwd(&dir_[0], bufferSize - 1);
#else
  if(getcwd(&dir_[0], bufferSize - 1) == nullptr) {
    TERM(-1);
  }
#endif
  dir_[bufferSize - 1] = '\0';
  const MString dir(&dir_[0]);

  // Get current execution command and command line arguments
  const MString executable(m_argv[0]);
  stringstream ss;
  if(m_argc > 1) {
    ss << MString(m_argv[0]);
    for(MInt n = 1; n < m_argc; n++) {
      ss << " " << MString(m_argv[n]);
    }
  } else {
    ss << "(none)";
  }
  const MString args(ss.str());

  // Create start timestamp
  array<char, bufferSize> dateString_{};
  time_t rawTime = 0;

  // Get the current time and write it to rawTime
  time(&rawTime);

  // Convert to time struct
  tm* timeInfo = localtime(&rawTime);

  // Format time to string and save to buffer
  strftime(&dateString_[0], bufferSize, "%Y-%m-%d %H:%M:%S", timeInfo);

  // Store date as real string
  const MString dateString(&dateString_[0]);

  // Print startup information
  cout << "      _____ ______    ________   ___   ________                             \n"
          " ____|\\   _ \\  _   \\ |\\   __  \\ |\\  \\ |\\   __  \\ ___                   \n"
          "  ___\\ \\  \\\\\\__\\ \\  \\\\ \\  \\|\\  \\\\ \\  \\\\ \\  \\|\\  \\ ___    \n"
          "   ___\\ \\  \\\\|__| \\  \\\\ \\   __  \\\\ \\  \\\\ \\   __  \\ ___         \n"
          "    ___\\ \\  \\    \\ \\  \\\\ \\  \\ \\  \\\\ \\  \\\\ \\  \\ \\  \\ ___    \n"
          "     ___\\ \\__\\    \\ \\__\\\\ \\__\\ \\__\\\\ \\__\\\\ \\__\\ \\__\\ ___    \n"
          "      ___\\|__|     \\|__| \\|__|\\|__| \\|__| \\|__|\\|__| ____               \n"
          "\n";
  cout << "Start time:            " << dateString << "\n"
       << "Number of ranks:       " << globalNoDomains() << "\n"
#ifdef _OPENMP
       << "Number of OMP threads: " << omp_get_max_threads() << "\n"
#endif
       << "Host (of rank 0):      " << host << "\n"
       << "Working directory:     " << dir << "\n"
       << "User:                  " << user << "\n"
       << "Executable:            " << executable << "\n"
       << "Command line args:     " << args << "\n"
       << endl;
}
