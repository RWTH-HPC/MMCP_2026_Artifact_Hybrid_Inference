// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_CONFIG_H_
#define MAIA_CONFIG_H_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file contains MAIA's configuration options.
///
/// Note: The following options can also be passed to the compiler with the
/// -D[option] (see also precompile script), e.g.:
/// \code
/// g++ -DMAIA_EXTRA_DEBUG ...
/// \endcode
///
/// To disable an option:
///   - comment it out here
///   - check out that precompile is not passing it to the make-file system
///   in the compilation mode that you are using.
////////////////////////////////////////////////////////////////////////////////
/// HARDWARE OPTIONS

/// HARDWARE PLATFORMS
// #define IBM_BLUE_GENE

/// \brief Used in Scratch to set alignment to cache line size, or not.
#define MAIA_SCRATCH_ALIGNMENT_BOUNDARY 64

////////////////////////////////////////////////////////////////////////////////
/// DEBUGGING AND PROFILING

/// \brief define to enable extra debugging checks (see precompile extra_debug)
/// #define MAIA_EXTRA_DEBUG

/// \brief enables asserts except for accessors
#if(defined(MAIA_EXTRA_DEBUG) || !defined(NDEBUG)) && !defined(MAIA_ASSERTS)
#define MAIA_ASSERTS
#endif

/// \brief enables asserts in accessors
#if defined(MAIA_EXTRA_DEBUG) && !defined(MAIA_ASSERT_ACCESSORS)
#define MAIA_ASSERT_ACCESSORS
#endif

/// \brief enables asserts in memory allocation
//#define MAIA_ASSERT_ALLOC

/// \brief Enable asserts in memory allocation in extra debug mode
#if defined(MAIA_EXTRA_DEBUG) && !defined(MAIA_ASSERT_ALLOC)
#define MAIA_ASSERT_ALLOC
#endif

/// \brief enables additional debug output in memory allocation
//#define MAIA_DEBUG_ALLOC

/// Enable allocation debug output in debug mode
#if !defined(NDEBUG) && !defined(MAIA_DEBUG_ALLOC)
#define MAIA_DEBUG_ALLOC
#endif

/// \brief Enable additional grid sanity checks (e.g. after balance/adaptation)
//#define MAIA_GRID_SANITY_CHECKS

/// Enable additional grid sanity checks in debug mode
#if !defined(NDEBUG) && !defined(MAIA_GRID_SANITY_CHECKS)
#define MAIA_GRID_SANITY_CHECKS
#endif

/// Enable debug output for DLB timer
// #define MAIA_DEBUG_DLB_TIMER

/// \brief define to enable profiling (incompatible with: MAIA_DEBUG_FUNCTION)
// #define MAIA_PROFILING

/// \brief define to enable debugging (incompatible with: MAIA_PROFILING)
// #define MAIA_DEBUG_FUNCTION

// \brief define to write massive parallel debug information to the m_log
//#define MP_DEBUG_LOG

// \brief define to enable use of access_properties file
//#define MAIA_WRITE_ACCESS_PROPERTIES_FILE

// \brief define to enable print out of properties; WARNING: will slow down the simulation
//#define MAIA_PRINT_PROPERTIES

/// \brief define to enable MAIA timers
#define MAIA_TIMER_FUNCTION

/// \brief Enable additional output for debugging timers
// #define MAIA_TIMER_DEBUG

/// \brief Assert that all timer checks pass, if not MAIA will terminate and print a debug trace
//#define MAIA_ASSERT_TIMER_CHECKS

/// \brief Enable additional checks when starting/stopping timers (enable for DEBUG, TIMER_DEBUG or
///        ASSERT_TIMER_CHECKS)
//#define MAIA_TIMER_CHECKS
#if(defined(MAIA_TIMER_DEBUG) || defined(MAIA_ASSERT_TIMER_CHECKS)) && !defined(MAIA_TIMER_CHECKS)
#define MAIA_TIMER_CHECKS
#endif

/// \brief synchronize timers across domains
//#define MAIA_SYNCHRONIZE_TIMERS

/// \brief define to use MPI_Wtime to get the time
#define MAIA_MPI_TIMER

/// \brief enables MPI debugging
/// #define MAIA_MPI_DEBUG

////////////////////////////////////////////////////////////////////////////////
/// LOGGING: InfoOut
/// For any questions regarding the usage, please refer to Michael (95188, mic@aia.rwth-aachen.de)
/// June 2012, mic

/// Set the project name you want in your m_log file in quotes, like "dummy_name".
#define MAIA_INFOOUT_PROJECT_NAME "n/a"

/* Choose the file type you want to use:
 * 0: Use a single file that is used by all processors, and that is written to using MPI I/O
 * 1: Use a physical file for each processor (as before, but with formatted output)
 *
 * This basically corresponds to the namespace MAIA_INFOOUT_FILETYPES in infoout.h
 */
#define MAIA_INFOOUT_FILE_TYPE 1

/* Set to 'true' if you want only the root domain (i.e. rank 0) to write to the log file,
 * otherwise to 'false'. This has different
 * meanings for the different file types:
 * 0: While all processors CAN write to the file, only rank 0 will do so.
 * 1: Only one file is created, and only rank 0 will write to it.
 */
#define MAIA_INFOOUT_ROOT_ONLY true

/* Set the minimum buffer size of the log file before output is written. Before
 * your log messages reach this size (in bytes), nothing is written to the file
 * (true for both file types). Especially for MPI files this can improve the
 * performance dramatically. Please note that of course the internal buffer will
 * be increased by at least this amount of bytes as well, resulting in a larger
 * memory footprint.
 *
 * If you set this to zero, each message (i.e. after each "endl") the buffer is
 * written to the file.
 * A good starting value is 256 kB (256*1024 = 262144).
 */
#define MAIA_INFOOUT_MIN_FLUSH_SIZE 0

////////////////////////////////////////////////////////////////////////////////
// The following options control optimizations that use compiler intrinsics

/// \brief define to enable (non-portable) intrinsic: __restrict
// #define USE_RESTRICT

/// \brief define to enables (non-portable) compiler attributes
// #define COMPILER_ATTRIBUTES

/// \brief define to enable (non-portable) assume aligned.
///
/// Note: this tells the compiler to _REALLY_ assume that "a" is aligned as a "b" even
/// if it cannot prove it. It returns immutable pointers by default, to get
/// mutable pointers use the _M variants.
// #define USE_ALIGNMENT

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// The following options enable/disable certain features from FV solver:

/// \brief disables FV Multigrid (makes the code slower and is not used)
// #define DISABLE_FV_MG

/// \brief enables Free Surface boundary condition:
// #define FV_FREE_SURFACE_BC

/// \brief enables marker which can be used to "mark" cells
/// it adds one additional int member to every cell for debugging purposes
/// and enables output of marker to solution file
//#define FV_MARKER

/// \brief This enables logging of access to the time-step variable when this
/// is not available yet.
#define MAIA_FV_LOG_ACCESS_TO_UNAVAILABLE_TIME_STEP

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// The following options enable/disable certain features in the DG solver

// Set how to handle the exchange of MPI surface data
//
// You need to specify one method of how to facilitate the exchange of MPI
// surface data. If 'DG_USE_MPI_BUFFERS' is enabled, data that needs to be
// exchanged is explicitly copied to contiguous buffers in memory for each
// domain. If 'DG_USE_MPI_DERIVED_TYPES' is enabled, individual MPI derived
// datatypes are created for each domain, and the (optional) buffering is left
// to the MPI implementation.
// Rationale: This switch is left in here until one method has been proven to
// be superior. Afterwards, the inferior method should be removed.
// Note: If neither or both methods are enabled, the code will not compiler.
#define DG_USE_MPI_BUFFERS
// #define DG_USE_MPI_DERIVED_TYPES

// If enabled, the DG solver will use MPI_Waitsome to mix receiving and memory
// copying in finishMpiExchange(). Will not compile if DG_USE_MPI_DERIVED_TYPES
// is enabled.
// #define DG_USE_MPI_WAITSOME
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// The following options control certain I/O aspects

/// Set Parallel NetCDF file type
/// (can be either NC_64BIT_OFFSET or NC_64BIT_DATA)
/// NC_64BIT_OFFSET: netcdf-3, 64 bit array length, matlab support
/// NC_64BIT_DATA: netcdf-5, 64 bit array length, 64 bit data types, no matlab support
#define MAIA_NCMPI_FILE_TYPE NC_64BIT_DATA

/// Enable fill mode for NetCDF variables, i.e. fill the variables with a default value (NaN/limits::max/min) first
/// before writing the actual data.
/// This will write all data twice to the file and increase IO overhead, however, it can be used to detect erroneous
/// files with incomplete written data (which has happened on Hawk) by checking for array fill values still present.
#define MAIA_NCMPI_FILL_VARIABLES false

/// Print NetCDF file hints after creating/opening and when closing files.
// #define MAIA_NCMPI_PRINT_FILE_HINTS

/// MPI-IO optimisation by setting striping etc.
#define MPI_IO_OPT

/// Print global MPI information at startup
#define MPI_IO_PRINT_INFO

/// Set the default backend class for ParallelIo if not defined on command line
/// Note: Under Windows, Parallel netCDF is not supported, thus in this case HDF5
///       is always chosen as the default backend.
#if defined(_WIN64)
#define PARALLELIO_DEFAULT_BACKEND ParallelIoHdf5
#endif
#if !defined(PARALLELIO_DEFAULT_BACKEND)
#define PARALLELIO_DEFAULT_BACKEND ParallelIoPNetcdf
#endif
////////////////////////////////////////////////////////////////////////////////

#endif
