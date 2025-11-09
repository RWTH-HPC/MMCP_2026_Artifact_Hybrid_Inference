// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "parallelio.h"
// needs to be included after parallelio.h!
#include "parallelio_pnetcdf.h"

#include "maiapnetcdf.h"

#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include "COMM/mpioverride.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "typetraits.h"

#if not defined(MAIA_MS_COMPILER)
#include <pwd.h>
#include <unistd.h>
#endif

#ifdef _SX
#include <sys/socket.h>
#include <sys/types.h>
#endif

using namespace maia;
using namespace parallel_io;
using namespace std;

// Check if valid Parallel NetCDF data file type was set
#if(MAIA_NCMPI_FILE_TYPE != NC_64BIT_OFFSET) && (MAIA_NCMPI_FILE_TYPE != NC_64BIT_DATA)
#error Bad value for MAIA_NCMPI_FILE_TYPE.
#endif

//------------------------------------------------------------------------------
// Parallel NetCDF-specific methods
//------------------------------------------------------------------------------

// Use unnamed namespace for some PnetCDF-specific magic (a.k.a. type traits)
namespace {

template <class DataType>
struct pnetcdf_traits {};

// MFloat traits
template <>
struct pnetcdf_traits<MFloat> {
  // Corresponding NetCDF data type
  static nc_type type() { return NC_DOUBLE; }

  // Write contiguously
  static int ncmpi_put_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MFloat* buf) {
    return ncmpi_put_vara_double_all(ncid, varid, start, count, buf);
  }

  // Write with stride
  static int ncmpi_put_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], const MFloat* buf) {
    return ncmpi_put_vars_double_all(ncid, varid, start, count, stride, buf);
  }

  // Read contiguously
  static int ncmpi_get_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     MFloat* buf) {
    return ncmpi_get_vara_double_all(ncid, varid, start, count, buf);
  }

  // Read with stride
  static int ncmpi_get_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], MFloat* buf) {
    return ncmpi_get_vars_double_all(ncid, varid, start, count, stride, buf);
  }

  // Set attribute
  static int ncmpi_put_att_type(int ncid, int varid, const char* name, nc_type xtype, const MPI_Offset nelems,
                                const MFloat* buf) {
    return ncmpi_put_att_double(ncid, varid, name, xtype, nelems, buf);
  }

  // Get attribute
  static int ncmpi_get_att_type(int ncid, int varid, const char* name, MFloat* buf) {
    return ncmpi_get_att_double(ncid, varid, name, buf);
  }
};

// MInt traits
template <>
struct pnetcdf_traits<MInt> {
  // Corresponding NetCDF data type
  static nc_type type() { return NC_INT; }

  // Write contiguously
  static int ncmpi_put_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MInt* buf) {
    return ncmpi_put_vara_int_all(ncid, varid, start, count, buf);
  }

  // Write with stride
  static int ncmpi_put_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], const MInt* buf) {
    return ncmpi_put_vars_int_all(ncid, varid, start, count, stride, buf);
  }

  // Read contiguously
  static int ncmpi_get_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     MInt* buf) {
    return ncmpi_get_vara_int_all(ncid, varid, start, count, buf);
  }

  // Read with stride
  static int ncmpi_get_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], MInt* buf) {
    return ncmpi_get_vars_int_all(ncid, varid, start, count, stride, buf);
  }

  // Set attribute
  static int ncmpi_put_att_type(int ncid, int varid, const char* name, nc_type xtype, MPI_Offset nelems,
                                const MInt* buf) {
    return ncmpi_put_att_int(ncid, varid, name, xtype, nelems, buf);
  }

  // Get attribute
  static int ncmpi_get_att_type(int ncid, int varid, const char* name, MInt* buf) {
    return ncmpi_get_att_int(ncid, varid, name, buf);
  }
};

// MLong traits
template <>
struct pnetcdf_traits<MLong> {
  // Corresponding NetCDF data type
  static nc_type type() { return NC_INT64; }

  // Write contiguously
  static int ncmpi_put_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MLong* buf) {
    return ncmpi_put_vara_long_all(ncid, varid, start, count, buf);
  }

  // Write with stride
  static int ncmpi_put_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], const MLong* buf) {
    return ncmpi_put_vars_long_all(ncid, varid, start, count, stride, buf);
  }

  // Read contiguously
  static int ncmpi_get_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     MLong* buf) {
    return ncmpi_get_vara_long_all(ncid, varid, start, count, buf);
  }

  // Read with stride
  static int ncmpi_get_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], MLong* buf) {
    return ncmpi_get_vars_long_all(ncid, varid, start, count, stride, buf);
  }

  // Set attribute
  static int ncmpi_put_att_type(int ncid, int varid, const char* name, nc_type xtype, MPI_Offset nelems,
                                const MLong* buf) {
    return ncmpi_put_att_long(ncid, varid, name, xtype, nelems, buf);
  }

  // Get attribute
  static int ncmpi_get_att_type(int ncid, int varid, const char* name, MLong* buf) {
    return ncmpi_get_att_long(ncid, varid, name, buf);
  }
};

// MChar traits
template <>
struct pnetcdf_traits<MChar> {
  // Corresponding NetCDF data type
  static nc_type type() { return NC_CHAR; }

  // Write contiguously
  static int ncmpi_put_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MChar* buf) {
    return ncmpi_put_vara_text_all(ncid, varid, start, count, buf);
  }

  // Write with stride
  static int ncmpi_put_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], const MChar* buf) {
    return ncmpi_put_vars_text_all(ncid, varid, start, count, stride, buf);
  }

  // Read contiguously
  static int ncmpi_get_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     MChar* buf) {
    return ncmpi_get_vara_text_all(ncid, varid, start, count, buf);
  }

  // Read with stride
  static int ncmpi_get_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], MChar* buf) {
    return ncmpi_get_vars_text_all(ncid, varid, start, count, stride, buf);
  }

  // Set attribute
  static int ncmpi_put_att_type(int ncid, int varid, const char* name, nc_type NotUsed(xtype), MPI_Offset nelems,
                                const MChar* buf) {
    // nc_type remains unused for _text API
    return ncmpi_put_att_text(ncid, varid, name, nelems, buf);
  }

  // Get attribute
  static int ncmpi_get_att_type(int ncid, int varid, const char* name, MChar* buf) {
    return ncmpi_get_att_text(ncid, varid, name, buf);
  }
};

// MUchar traits
template <>
struct pnetcdf_traits<MUchar> {
  // Corresponding NetCDF data type
  static nc_type type() { return NC_UBYTE; }

  // Write contiguously
  static int ncmpi_put_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MUchar* buf) {
    return ncmpi_put_vara_uchar_all(ncid, varid, start, count, buf);
  }

  // Write with stride
  static int ncmpi_put_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], const MUchar* buf) {
    return ncmpi_put_vars_uchar_all(ncid, varid, start, count, stride, buf);
  }

  // Read contiguously
  static int ncmpi_get_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     MUchar* buf) {
    return ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf);
  }

  // Read with stride
  static int ncmpi_get_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], MUchar* buf) {
    return ncmpi_get_vars_uchar_all(ncid, varid, start, count, stride, buf);
  }

  // Set attribute
  static int ncmpi_put_att_type(int ncid, int varid, const char* name, nc_type xtype, MPI_Offset nelems,
                                const MUchar* buf) {
    return ncmpi_put_att_uchar(ncid, varid, name, xtype, nelems, buf);
  }

  // Get attribute
  static int ncmpi_get_att_type(int ncid, int varid, const char* name, MUchar* buf) {
    return ncmpi_get_att_uchar(ncid, varid, name, buf);
  }
};

// MUlong traits
template <>
struct pnetcdf_traits<MUlong> {
  // Corresponding NetCDF data type
  static nc_type type() { return NC_UINT64; }

  // Write contiguously
  static int ncmpi_put_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MUlong* buf) {
    return ncmpi_put_vara_ulonglong_all(ncid, varid, start, count, (const long long unsigned int*)(buf));
  }

  // Write with stride
  static int ncmpi_put_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], const MUlong* buf) {
    return ncmpi_put_vars_ulonglong_all(ncid, varid, start, count, stride, (const long long unsigned int*)(buf));
  }

  // Read contiguously
  static int ncmpi_get_vara_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     MUlong* buf) {
    return ncmpi_get_vara_ulonglong_all(ncid, varid, start, count, (long long unsigned int*)(buf));
  }

  // Read with stride
  static int ncmpi_get_vars_type_all(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                     const MPI_Offset stride[], MUlong* buf) {
    return ncmpi_get_vars_ulonglong_all(ncid, varid, start, count, stride, (long long unsigned int*)(buf));
  }

  // Set attribute
  static int ncmpi_put_att_type(int ncid, int varid, const char* name, nc_type xtype, MPI_Offset nelems,
                                const MUlong* buf) {
    return ncmpi_put_att_ulonglong(ncid, varid, name, xtype, nelems, (const long long unsigned int*)(buf));
  }

  // Get attribute
  static int ncmpi_get_att_type(int ncid, int varid, const char* name, MUlong* buf) {
    return ncmpi_get_att_ulonglong(ncid, varid, name, (long long unsigned int*)(buf));
  }
};

} // namespace

//
//------------------------------------------------------------------------------
// Static file system-related methods
//------------------------------------------------------------------------------
/**
 * \brief Check if specified file is a valid HDF5 file (i.e. can be
 *        opened). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name Name of the file that is to be checked.
 * \param[in] mpiComm MPI communicator to be used.
 *
 * \details This must be a collective call from all ranks in mpiComm.
 */
MBool ParallelIoPNetcdf::b_isValidFile(const MString& name, const MPI_Comm& mpiComm) {
  TRACE();

  MBool returnValue;

  MInt status, fileId;
  status = ncmpi_open(mpiComm, name.c_str(), NC_NOWRITE, globalMpiInfo(), &fileId);

  if(status == NC_NOERR) {
    returnValue = true;
    status = ncmpi_close(fileId);
    b_error(status, name, AT_);
  } else {
    returnValue = false;
  }

  return returnValue;
}


/**
 * \brief Returns backend-specific ending of filename (either ".Netcdf" or
 *  ".Hdf5")
 * \author Konstantin Froehlich
 * \date 2015-10-21
 *
 */
const MString ParallelIoPNetcdf::b_fileExt() {
  TRACE();

  return ".Netcdf";
}


//------------------------------------------------------------------------------
// Constructor & Destructor
//------------------------------------------------------------------------------
/**
 * \brief Creates a new object to read and write *big* data files. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-2
 * \author (Modified by: )Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>
 * \date 2015-03-01
 *
 * \param[in] fileMode The file mode can be either maia::parallel_io::PIO_CREATE,
 *                     maia::parallel_io::PIO_APPEND, or maia::parallel_io::PIO_READ.
 * \param[in] mpiComm The MPI communicator that should be used to open/create
 *                    the file.
 */
ParallelIoPNetcdf::ParallelIoPNetcdf(const MString& fileName, MInt fileMode, const MPI_Comm& mpiComm)
  : ParallelIoBase<ParallelIoPNetcdf>(fileName, fileMode, mpiComm) {
  TRACE();

#ifdef DISABLE_OUTPUT
  if(m_fileMode != PIO_READ) return;
#endif

  switch(m_fileMode) {
    case PIO_CREATE: {
      // Create a new file (do not overwrite existing)
      MInt status =
          ncmpi_create(m_mpiComm, m_fileName.c_str(), NC_NOCLOBBER | MAIA_NCMPI_FILE_TYPE, globalMpiInfo(), &b_ncId);
      b_error(status, m_fileName, AT_);
    } break;

    case PIO_REPLACE: {
      // Create a new file (overwrite existing)
      MInt status =
          ncmpi_create(m_mpiComm, m_fileName.c_str(), NC_CLOBBER | MAIA_NCMPI_FILE_TYPE, globalMpiInfo(), &b_ncId);
      b_error(status, m_fileName, AT_);
    } break;

    case PIO_APPEND: {
      // Attempt to open an existing file to append data
      MInt status =
          ncmpi_open(m_mpiComm, m_fileName.c_str(), NC_WRITE | MAIA_NCMPI_FILE_TYPE, globalMpiInfo(), &b_ncId);
      b_error(status, m_fileName, AT_);

      // Set correct data mode status
      b_ncDataMode = true;

      // Read in all existing dimensions and populate dimensions map
      int noDims;
      status = ncmpi_inq_ndims(b_ncId, &noDims);
      b_error(status, m_fileName, AT_);
      for(int dimId = 0; dimId < noDims; dimId++) {
        MPI_Offset dimLength;
        status = ncmpi_inq_dimlen(b_ncId, dimId, &dimLength);
        b_error(status, m_fileName, AT_);
        if(b_ncDimensions.count(dimLength) == 0u) {
          b_ncDimensions[dimLength] = NcDimProxy();
          b_ncDimensions[dimLength].dimNo = dimLength;
          b_ncDimensions[dimLength].dimId = dimId;
        }
      }

      // Enter definition mode so that new attributes/variables may be added
      b_ncRedef();
    } break;

    case PIO_READ: {
      // Open file for reading
      MInt status = ncmpi_open(m_mpiComm, m_fileName.c_str(), NC_NOWRITE, globalMpiInfo(), &b_ncId);
      b_error(status, m_fileName, AT_);

      // Set correct data mode status
      b_ncDataMode = true;
    } break;

    default: {
      mTerm(1, AT_, "Unsupported file mode.");
    } break;
  }

#ifdef MAIA_NCMPI_PRINT_FILE_HINTS
  if(m_domainId == 0) {
    std::cerr << std::endl << "Created/replaced/opened file: " << m_fileName << std::endl;
    b_printFileHints();
  }
#endif
}


/**
 * \brief Calls close(). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-26
 */
ParallelIoPNetcdf::~ParallelIoPNetcdf() {
  TRACE();

  // Check if file is already closed
  if(b_ncId != -1) {
    close();
  }
}


/// \brief Close the file (normally called by the destructor but needs to be explicitely called earlier in specific
/// cases, e.g. when
// the file MPI communicator is destroyed before the file object goes out of scope).
void ParallelIoPNetcdf::close() {
  TRACE();
#ifdef MAIA_NCMPI_PRINT_FILE_HINTS
  if(m_domainId == 0) {
    std::cerr << std::endl << "Closing file: " << m_fileName << std::endl;
    b_printFileHints();
  }
#endif

  MInt status = ncmpi_close(b_ncId);
  b_error(status, m_fileName, AT_);
  b_ncId = -1; // Reset Netcdf file id, calling any other function after close() should result in an invalid file id
               // Netcdf error

#ifdef MAIA_EXTRA_DEBUG
  for(set<MString>::iterator it = m_unwrittenArrays.begin(); it != m_unwrittenArrays.end(); ++it) {
    cerr << "Warning: array '" << *it << "' in file '" << m_fileName << "' was defined but never written. "
         << "Make sure that this is the intended behavior." << endl;
  }
  for(set<MString>::iterator it = m_unwrittenScalars.begin(); it != m_unwrittenScalars.end(); ++it) {
    cerr << "Warning: scalar '" << *it << "' in file '" << m_fileName << "' was defined but never written. "
         << "Make sure that this is the intended behavior." << endl;
  }
#endif
}

//------------------------------------------------------------------------------
// File methods
//------------------------------------------------------------------------------


/**
 * \brief Leave the define mode (NetCDF-specific). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-26
 */
void ParallelIoPNetcdf::b_ncEndDef() {
  TRACE();

  if(!b_ncDataMode) {
    MInt status = ncmpi_enddef(b_ncId);
    b_error(status, m_fileName, AT_);
    b_ncDataMode = true;
  }
}


/**
 * \brief Enter define mode (NetCDF-specific). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-26
 */
void ParallelIoPNetcdf::b_ncRedef() {
  TRACE();

  if(b_ncDataMode) {
    MInt status = ncmpi_redef(b_ncId);
    b_error(status, m_fileName, AT_);
    b_ncDataMode = false;
  }
}


/**
 * \brief Adds all additional header information that are needed in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-09
 *
 * \details This method must be called exactly once for newly created files.
 */
void ParallelIoPNetcdf::b_saveHeader() {
  TRACE();

  b_addAdditionalHeader();
  b_writeAdditionalData();
}


// Forward declaration of specialization for use in b_addAdditionalHeader
template <>
void ParallelIoPNetcdf::b_setAttribute(const MString* value,
                                       const MString& name,
                                       const MString& datasetName,
                                       const size_type totalCount);


/**
 * \brief Write additional headers to file (e.g. grid file name, creation date
 *        etc.). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-12
 */
void ParallelIoPNetcdf::b_addAdditionalHeader() {
  TRACE();

  b_ncRedef();

  // For a newly created or modified file, add some meta information
  if(m_fileMode != PIO_READ) {
    const MInt maxNoChars = 256;

    // Get all meta-data on root process & communicate, since some of the C
    // functions might not be thread-safe (at least getpwuid() is not)
    MChar user[maxNoChars];
    MChar host[maxNoChars];
    MChar dir[maxNoChars];
    MChar exec[maxNoChars];
    MChar date[maxNoChars];

    // Create object & initialize data
    fill(user, user + maxNoChars, '\0');
    fill(host, host + maxNoChars, '\0');
    fill(dir, dir + maxNoChars, '\0');
    fill(exec, exec + maxNoChars, '\0');
    fill(date, date + maxNoChars, '\0');

    if(m_domainId == 0) {
      // Gets the current username
      passwd* p;
      p = getpwuid(getuid());
      if(p) {
        strncpy(user, p->pw_name, maxNoChars - 1);
      } else {
        strncpy(user, "n/a", maxNoChars - 1);
      }

      // Gets the current hostname
      gethostname(host, maxNoChars - 1);
      strcpy(&host[strlen(host)], " (");
      strcpy(&host[strlen(host)], XSTRINGIFY(MAIA_HOST_STRING));
      strcpy(&host[strlen(host)], ")");

// Gets the current directory
#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#endif

      getcwd(dir, maxNoChars - 1);

      readlink("/proc/self/exe", exec, maxNoChars - 1);

#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic pop
#endif

      // Get the current time and write it to rawTime
      time_t rawTime;
      time(&rawTime);

      // Convert to time struct
      tm* timeInfo;
      timeInfo = localtime(&rawTime);

      // Format time to string and save to buffer
      strftime(date, maxNoChars, "%Y-%m-%d %H:%M:%S", timeInfo);
    }

    // Pack data
    const MInt noItems = 5;
    MChar buffer[noItems * maxNoChars];
    memcpy(buffer + 0 * maxNoChars, user, maxNoChars);
    memcpy(buffer + 1 * maxNoChars, host, maxNoChars);
    memcpy(buffer + 2 * maxNoChars, dir, maxNoChars);
    memcpy(buffer + 3 * maxNoChars, exec, maxNoChars);
    memcpy(buffer + 4 * maxNoChars, date, maxNoChars);

    // Broadcast time from rank 0 to ensure that every rank has the same
    // information
    MPI_Bcast(&buffer, noItems * maxNoChars, MPI_CHAR, 0, m_mpiComm, AT_, "buffer");

    // Unpack data
    memcpy(user, buffer + 0 * maxNoChars, maxNoChars);
    memcpy(host, buffer + 1 * maxNoChars, maxNoChars);
    memcpy(dir, buffer + 2 * maxNoChars, maxNoChars);
    memcpy(exec, buffer + 3 * maxNoChars, maxNoChars);
    memcpy(date, buffer + 4 * maxNoChars, maxNoChars);
    MString version = MString(XSTRINGIFY(MAIA_VERSION_STRING));
    MString build = MString(XSTRINGIFY(MAIA_COMPILER_STRING)) + " " + MString(XSTRINGIFY(MAIA_BUILD_TYPE_STRING)) + " ("
                    + MString(XSTRINGIFY(MAIA_COMPILER_VERSION_STRING)) + ")";

    if(m_fileMode == PIO_CREATE || m_fileMode == PIO_REPLACE) {
      // Add file attributes only needed for creation
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(user), "_meta_creation_user");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(host), "_meta_creation_host");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(dir), "_meta_creation_directory");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(exec), "_meta_creation_executable");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(m_noDomains, "_meta_creation_noDomains");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(date), "_meta_creation_date");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(version, "_meta_creation_revision");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(build, "_meta_creation_build");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(ncmpi_inq_libvers(), "_meta_creation_pnetcdf_version");
    } else if(m_fileMode == PIO_APPEND) {
      // Add file attributes that should be set at each modification
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(user), "_meta_lastModified_user");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(host), "_meta_lastModified_host");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(dir), "_meta_lastModified_directory");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(exec), "_meta_lastModified_executable");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(m_noDomains, "_meta_lastModified_noDomains");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(MString(date), "_meta_lastModified_date");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(version, "_meta_lastModified_revision");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(build, "_meta_lastModified_build");
      ParallelIoBase<ParallelIoPNetcdf>::setAttribute(ncmpi_inq_libvers(), "_meta_lastModified_pnetcdf_version");
    }
  }
}


/**
 * \brief Write additional data to file (NetCDF-specific). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-26
 */
void ParallelIoPNetcdf::b_writeAdditionalData() {
  TRACE();

  // At the moment, nothing happens here
}

void ParallelIoPNetcdf::b_defineArray(maiabd_type type, const MString& name, const MString& path,
                                      const size_type noDims, const size_type* totalCount) {
  TRACE();
  (void)type;
  (void)name;
  (void)path;
  (void)noDims;
  (void)totalCount;
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
}


//------------------------------------------------------------------------------
// Define mode methods
//------------------------------------------------------------------------------
/**
 * \brief Defines an array in the file. <b>[MPI]</b>
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar
 * \param[in] type Data type of the array (may be either
 *                 maia::parallel_io::PIO_FLOAT or maia::parallel_io::PIO_INT or
 *                 maia::parallel_io::PIO_STRING or maia::parallel_io::PIO_UCHAR).
 * \param[in] name Name of the array (must not be empty).
 * \param[in] noDims Number of array dimensions
 * \param[in] totalCount Total size of the array in each dimension.
 *
 */
void ParallelIoPNetcdf::b_defineArray(maiabd_type type, const MString& name, size_type noDims, size_type* totalCount) {
  TRACE();

  b_ncRedef();

  // Fill values
  MFloat fillValueFloat = MFloatNaN;
  MInt fillValueInt = std::numeric_limits<MInt>::max();
  MLong fillValueLong = std::numeric_limits<MLong>::max();
  MUlong fillValueUlong = std::numeric_limits<MUlong>::max();
  MUchar fillValueChar = std::numeric_limits<MChar>::max();
  MUchar fillValueUchar = std::numeric_limits<MUchar>::max();

  // Determine NC data type
  nc_type dataType;
  [[maybe_unused]] void* fillValue = nullptr;
  if(type == PIO_FLOAT) {
    dataType = NC_DOUBLE;
    fillValue = &fillValueFloat;
  } else if(type == PIO_INT) {
    dataType = NC_INT;
    fillValue = &fillValueInt;
  } else if(type == PIO_LONG) {
    dataType = NC_INT64;
    fillValue = &fillValueLong;
  } else if(type == PIO_STRING) {
    dataType = NC_CHAR;
    fillValue = &fillValueChar;
  } else if(type == PIO_UCHAR) {
    dataType = NC_UBYTE;
    fillValue = &fillValueUchar;
  } else if(type == PIO_ULONGLONG) {
    dataType = NC_UINT64;
    fillValue = &fillValueUlong;
  } else {
    TERMM(1, "Invalid ParallelIo data type!");
  }

  // Determine whether (a) new dimension(s) needs to be created and if yes,
  // create it
  MInt dimId, status, varId;
  MIntScratchSpace dimIds(noDims, FUN_, "dimIds");
  for(size_type dId = 0; dId < noDims; dId++) {
    if(b_ncDimensions.count(totalCount[dId]) != 0u) {
      // Dimension was found, now get its id
      dimId = b_ncDimensions.find(totalCount[dId])->second.dimId;
    } else {
      // Get next dimension name
      MInt maxUsedDimensionNo = -1;
      for(NcDimMap::const_iterator it = b_ncDimensions.begin(); it != b_ncDimensions.end(); it++) {
        maxUsedDimensionNo = max(maxUsedDimensionNo, it->second.dimNo);
      }

      // Set dimension name according to found dimension numbers
      const MInt dimNo = maxUsedDimensionNo + 1;
      const MString dimensionName = "dim" + to_string(dimNo);

      // Define new dimension in the file and store its information
      status = ncmpi_def_dim(b_ncId, dimensionName.c_str(), totalCount[dId], &dimId);
      b_error(status, dimensionName, AT_);

      // Add dimension to map
      b_ncDimensions[totalCount[dId]] = NcDimProxy();
      b_ncDimensions[totalCount[dId]].dimNo = dimNo;
      b_ncDimensions[totalCount[dId]].dimId = dimId;
    }
    dimIds[dId] = dimId;
  }

  // Define dimension and variable in the file
  status = ncmpi_def_var(b_ncId, name.c_str(), dataType, noDims, &dimIds[0], &varId);
  b_error(status, name, AT_);

#if MAIA_NCMPI_FILL_VARIABLES
  // NOTE: set fill mode and default fill value to detect unwritten variables/detect IO errors
  status = ncmpi_def_var_fill(b_ncId, varId, 0, fillValue);
  b_error(status, name, AT_);
#endif
}


/**
 * \brief Defines a scalar in the file. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] type Data type of the scalar (may be either
 *                 maia::parallel_io::PIO_FLOAT or maia::parallel_io::PIO_INT).
 * \param[in] name Name of the scalar (must not be empty).
 */
void ParallelIoPNetcdf::b_defineScalar(maiabd_type type, const MString& name) {
  TRACE();

  // Choose the HDF5 specific native datatype or a string of variable length
  b_ncRedef();

  // Fill values
  MFloat fillValueFloat = MFloatNaN;
  MInt fillValueInt = std::numeric_limits<MInt>::max();
  MLong fillValueLong = std::numeric_limits<MLong>::max();
  MUlong fillValueUlong = std::numeric_limits<MUlong>::max();
  MUchar fillValueChar = std::numeric_limits<MChar>::max();
  MUchar fillValueUchar = std::numeric_limits<MUchar>::max();

  // Determine NC data type
  nc_type dataType;
  [[maybe_unused]] void* fillValue = nullptr;
  if(type == PIO_FLOAT) {
    dataType = NC_DOUBLE;
    fillValue = &fillValueFloat;
  } else if(type == PIO_INT) {
    dataType = NC_INT;
    fillValue = &fillValueInt;
  } else if(type == PIO_LONG) {
    dataType = NC_INT64;
    fillValue = &fillValueLong;
  } else if(type == PIO_STRING) {
    dataType = NC_CHAR;
    fillValue = &fillValueChar;
  } else if(type == PIO_UCHAR) {
    dataType = NC_UBYTE;
    fillValue = &fillValueUchar;
  } else if(type == PIO_ULONGLONG) {
    dataType = NC_UINT64;
    fillValue = &fillValueUlong;
  } else {
    TERMM(1, "Invalid ParallelIo data type!");
  }

  // For a scalar, only the variable needs to be defined
  MInt varId;
  MInt status = ncmpi_def_var(b_ncId, name.c_str(), dataType, 0, nullptr, &varId);
  b_error(status, name, AT_);

#if MAIA_NCMPI_FILL_VARIABLES
  // NOTE: set fill mode and default fill value to detect unwritten variables/detect IO errors
  status = ncmpi_def_var_fill(b_ncId, varId, 0, fillValue);
  b_error(status, name, AT_);
#endif
}


//------------------------------------------------------------------------------
// Inquiry methods
//------------------------------------------------------------------------------
/**
 * \brief Check if the file contains a dataset with the given name and
 *        number of dimensions.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name The name of the dataset that should be checked.
 * \param[in] noDimensions Number of dimensions of the dataset to match,
 *                         default -1 which is any dimensionality, 0 will only
 *                         check scalars, 1 arrays with one dimension, ...
 *
 * \return True if an dataset of the given name with the given number of
 *         dimensions exists.
 */
MBool ParallelIoPNetcdf::b_hasDataset(const MString& name, const size_type noDimensions) {
  TRACE();

  // Get number of variable in file
  MInt noVars;
  MInt status = ncmpi_inq_nvars(b_ncId, &noVars);
  b_error(status, m_fileName, AT_);

  // Check each variable if
  // - names match
  // - has the given dimension (or match any dimension if -1 which is default)
  MBool varExists = false;
  for(MInt i = 0; i < noVars; i++) {
    MChar varname[NC_MAX_NAME + 1];
    status = ncmpi_inq_varname(b_ncId, i, varname);
    b_error(status, m_fileName, AT_);

    MInt nDims;
    status = ncmpi_inq_varndims(b_ncId, i, &nDims);
    b_error(status, varname, AT_);

    if(name == varname && (nDims == noDimensions || noDimensions == -1)) {
      varExists = true;
      break;
    }
  }

  return varExists;
}

MBool ParallelIoPNetcdf::b_hasObject(const MString& path) {
  TRACE();

  // pnetcdf doesn't support groups (hdf5 feature), therefore this
  // function can only check for datasets
  return b_hasDataset(path, -1);
}

MBool ParallelIoPNetcdf::b_hasDataset(const MString& name, const MString& path) {
  TRACE();

  (void)path;
  (void)name;
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
  return false;
}

void ParallelIoPNetcdf::b_getDatasetNames(std::vector<MString>& names, const MString& path) {
  TRACE();

  (void)path;
  (void)names;
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
}

void ParallelIoPNetcdf::b_getGroupNames(std::vector<MString>& names, const MString& path) {
  TRACE();

  (void)path;
  (void)names;
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
}

ParallelIo::size_type ParallelIoPNetcdf::b_getDatasetNoDims(const MString& path, const MString& name) {
  TRACE();

  (void)path;
  (void)name;
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
  return -1;
}

void ParallelIoPNetcdf::b_getDatasetSize(const MString& name, const MString& path, size_type noDims, size_type* data) {
  TRACE();

  (void)path;
  (void)name;
  (void)noDims;
  (void)data;
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
}


/**
 * \brief Returns the data type of a dataset in the file (can be array, multi-D
 * array or scalar).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name Name of the dataset to check.
 *
 * \return The data type as defined in the namespace maia::parallel_io.
 */
MInt ParallelIoPNetcdf::b_getDatasetType(const MString& name) {
  TRACE();

  // Get variable id
  MInt varId;
  MInt status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Get variable type
  nc_type ncType;
  status = ncmpi_inq_vartype(b_ncId, varId, &ncType);
  b_error(status, name, AT_);

  // Translate NC type to ParallelIo type
  MInt typeId;
  if(ncType == NC_INT) {
    typeId = PIO_INT;
  } else if(ncType == NC_DOUBLE) {
    typeId = PIO_FLOAT;
  } else if(ncType == NC_INT64) {
    typeId = PIO_LONG;
  } else if(ncType == NC_CHAR) {
    typeId = PIO_STRING;
  } else if(ncType == NC_UBYTE) {
    typeId = PIO_UCHAR;
  } else if(ncType == NC_UINT64) {
    typeId = PIO_ULONGLONG;
  } else {
    typeId = PIO_UNKNOWN_TYPE;
  }

  return typeId;
}


/**
 * \brief Returns a vector of all existing datasets with the given number of
 * dimensions in the file (if any).
 *        <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] dimension Match only datasets with given number of dimensions,
 *                      default -1 which will return all datasets, 0 will only
 *                      return scalars, 1 arrays with one dimension, ...
 * \param[out] names List of dataset names.
 */
void ParallelIoPNetcdf::b_getDatasetNames(vector<MString>& names, const size_type dimension) {
  TRACE();

  // Erase vector contents
  vector<MString>().swap(names);

  MInt status;

  // Get number of variable in file
  MInt noVars;
  status = ncmpi_inq_nvars(b_ncId, &noVars);
  b_error(status, m_fileName, AT_);

  // Check each variable if
  // - has the specified number of dimensions or return all if dimension=-1
  // (default)
  for(MInt i = 0; i < noVars; i++) {
    MInt noDimensions;
    status = ncmpi_inq_varndims(b_ncId, i, &noDimensions);
    b_error(status, m_fileName, AT_);

    if(noDimensions == dimension || dimension == -1) {
      // If this is a data file, get data file-specific name, otherwise just get
      // the variable name
      MChar varname[NC_MAX_NAME + 1];
      status = ncmpi_inq_varname(b_ncId, i, varname);
      b_error(status, m_fileName, AT_);
      names.emplace_back(varname);
    }
  }
}


/**
 * \brief Get number of dimensions of a dataset with the given name.
 *
 * \author Konstantin Froehlich
 * \date 2015-10-09
 *
 * \param[in] name Dataset name for which the number of dimensions is returned
 *
 * The number of dimensions of a dataset is returned, 0 is a scalar, 1 an array
 * with one dimension, ...
 *
 * \return Number of dimensions of the specified dataset
 */
ParallelIo::size_type ParallelIoPNetcdf::b_getDatasetNoDims(const MString& name) {
  TRACE();

  // Get variable id
  MInt varId;
  MInt status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Get number of variable dimensions
  MInt noDims;
  status = ncmpi_inq_varndims(b_ncId, varId, &noDims);
  b_error(status, name, AT_);

  return noDims;
}


/**
 * \brief Get the length of one dimension of an arbitrary array in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name Name of the array to check.
 * \param[in] dimensionId dimension of the array to check.
 *
 * \return Total number of elements of the array in the given dimension.
 */
ParallelIo::size_type ParallelIoPNetcdf::b_getDatasetSize(const MString& name, const size_type dimensionId) {
  TRACE();

  // Get variable id
  MInt varId;
  MInt status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Get number of array dimensions
  size_type noDims = b_getDatasetNoDims(name);

  // Get variable dimension
  MIntScratchSpace dimId(noDims, FUN_, "dimId");
  status = ncmpi_inq_vardimid(b_ncId, varId, &dimId[0]);
  b_error(status, name, AT_);

  // Get variable size
  MPI_Offset arraySize;
  status = ncmpi_inq_dimlen(b_ncId, dimId[dimensionId], &arraySize);
  b_error(status, name, AT_);

  return static_cast<size_type>(arraySize);
}


/**
 * \brief Check if a given attribute exists in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name Name of the attribute to be checked.
 * \param[in] datasetName Name of the dataset for which the attribute is
 *                        checked. If empty, a file attribute is checked.
 *
 * \return True, if the attribute exists.
 */
MBool ParallelIoPNetcdf::b_hasAttribute(const MString& name, const MString& datasetName) {
  TRACE();

  // Determine variable id
  MInt varId;
  if(datasetName.empty()) {
    varId = NC_GLOBAL;
  } else {
    MInt status = ncmpi_inq_varid(b_ncId, datasetName.c_str(), &varId);
    b_error(status, datasetName, AT_);
  }

  // Determine if attribute exists
  MInt attId;
  MInt status = ncmpi_inq_attid(b_ncId, varId, name.c_str(), &attId);

  // Generate return falue
  MBool attributeExists;
  if(status == NC_NOERR) {
    attributeExists = true;
  } else if(status == NC_ENOTATT) {
    attributeExists = false;
  } else {
    attributeExists = false;
    b_error(status, name, AT_);
  }

  return attributeExists;
}


/**
 * \brief Returns the data type of an attribute in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name Name of the attribute to be checked.
 * \param[in] datasetName Name of the dataset for which the attribute is
 *                        checked. If empty, a file attribute is checked.
 *
 * \return The data type as defined in the namespace maia::parallel_io.
 */
MInt ParallelIoPNetcdf::b_getAttributeType(const MString& name, const MString& datasetName) {
  TRACE();

  // Determine variable id
  MInt varId;
  if(datasetName.empty()) {
    varId = NC_GLOBAL;
  } else {
    MInt status = ncmpi_inq_varid(b_ncId, datasetName.c_str(), &varId);
    b_error(status, datasetName, AT_);
  }

  // Get attribute type
  nc_type ncType;
  MInt status = ncmpi_inq_atttype(b_ncId, varId, name.c_str(), &ncType);
  b_error(status, name, AT_);

  // Translate NC type to ParallelIo type
  MInt typeId;
  if(ncType == NC_INT) {
    typeId = PIO_INT;
  } else if(ncType == NC_DOUBLE) {
    typeId = PIO_FLOAT;
  } else if(ncType == NC_INT64) {
    typeId = PIO_LONG;
  } else if(ncType == NC_CHAR) {
    typeId = PIO_STRING;
  } else if(ncType == NC_UBYTE) {
    typeId = PIO_UCHAR;
  } else if(ncType == NC_UINT64) {
    typeId = PIO_ULONGLONG;
  } else {
    typeId = PIO_UNKNOWN_TYPE;
  }

  return typeId;
}


/**
 * \brief Print PNetCDF file hints to cerr.
 *
 * Source: https://trac.mcs.anl.gov/projects/parallel-netcdf/browser/trunk/examples/C/hints.c
 */
void ParallelIoPNetcdf::b_printFileHints() {
  TRACE();

  char value[MPI_MAX_INFO_VAL];
  MInt status, len, flag;
  MPI_Offset header_size, header_extent;
  MPI_Offset h_align = -1, v_align = -1, h_chunk = -1;
  MPI_Info info_used;

  // Get header size
  status = ncmpi_inq_header_size(b_ncId, &header_size);
  b_error(status, m_fileName, AT_);

  // Get maximum size of header (without the need to move data)
  status = ncmpi_inq_header_extent(b_ncId, &header_extent);
  b_error(status, m_fileName, AT_);

  // Get file MPI information
  status = ncmpi_inq_file_info(b_ncId, &info_used);
  b_error(status, m_fileName, AT_);

  // Get header align size (in bytes)
  MPI_Info_get_valuelen(info_used, "nc_header_align_size", &len, &flag, AT_);
  if(flag != 0) {
    MPI_Info_get(info_used, "nc_header_align_size", len + 1, value, &flag, AT_);
    h_align = strtoll(value, nullptr, 10);
  }

  // Get variable align size (in bytes)
  MPI_Info_get_valuelen(info_used, "nc_var_align_size", &len, &flag, AT_);
  if(flag != 0) {
    MPI_Info_get(info_used, "nc_var_align_size", len + 1, value, &flag, AT_);
    v_align = strtoll(value, nullptr, 10);
  }

  // Get header read chuck size (in bytes)
  MPI_Info_get_valuelen(info_used, "nc_header_read_chunk_size", &len, &flag, AT_);
  if(flag != 0) {
    MPI_Info_get(info_used, "nc_header_read_chunk_size", len + 1, value, &flag, AT_);
    h_chunk = strtoll(value, nullptr, 10);
  }

  MPI_Info_free(&info_used, AT_);

  // Output file hint information
  std::cerr << "##### PNetCDF file hints #####" << std::endl;

  if(h_align == -1) {
    std::cerr << "nc_header_align_size      is NOT set" << std::endl;
  } else {
    std::cerr << "nc_header_align_size      set to = " << h_align << std::endl;
  }

  if(v_align == -1) {
    std::cerr << "nc_var_align_size         is NOT set" << std::endl;
  } else {
    std::cerr << "nc_var_align_size         set to = " << v_align << std::endl;
  }
  if(h_chunk == -1) {
    std::cerr << "nc_header_read_chunk_size is NOT set" << std::endl;
  } else {
    std::cerr << "nc_header_read_chunk_size set to = " << h_chunk << std::endl;
  }

  std::cerr << "header size                      = " << header_size << std::endl;
  std::cerr << "header extent                    = " << header_extent << std::endl;

  // Check free header space in append mode
  if(m_fileMode == PIO_APPEND) {
    const MPI_Offset header_free = header_extent - header_size;
    std::cerr << "header free space (append mode)  = " << header_free << std::endl;
    if(header_free < 1024) {
      std::cerr << "WARNING: ParallelIoPNetcdf file in append mode has less than 1KB of free header "
                   "space. If this space is used up by adding new data to the header the entire data "
                   "file has to be moved to make room for the new definitions which may cause MPI I/O "
                   "errors."
                << std::endl;
    }
  }

  std::cerr << "##### PNetCDF file hints #####" << std::endl;
}


template <class T>
void ParallelIoPNetcdf::b_writeArray(const T* array, const MString& path, const MString& name, const size_type noDims,
                                     const size_type* start, const size_type* count, const size_type* ghost) {
  (void)array;
  (void)path;
  (void)name;
  (void)noDims;
  (void)start;
  (void)count;
  (void)ghost;
  TRACE();
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
}


//------------------------------------------------------------------------------
// Data mode methods
//------------------------------------------------------------------------------
/**
 * \brief Writes array data to file (generic version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the data to write.
 * \param[in] array Pointer to data in memory.
 * \param[in] name Name of the dataset.
 * \param[in] noDims Number of array dimensions
 * \param[in] start Start offset in the *file* (i.e. the domain should start
 *                  writing from here).
 * \param[in] count Number of elements written from this domain.
 * \param[in] memoryStride Stride in *memory* between data values.
 */
template <class T>
void ParallelIoPNetcdf::b_writeArray(const T* const array, const MString& name, const size_type noDims,
                                     MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                     const size_type noChunks, MPI_Offset diskStride) {
  TRACE();

  b_ncEndDef();

  // Get variable id
  MInt varId;
  MInt status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine total data count
  size_type totalCount = 1;
  for(size_type d = 0; d < noDims; d++) {
    totalCount *= count[d];
  }

  // Create temporary storage space if needed and set data pointers
  MInt tmpScratchSize = (memoryStride == 1) ? 1 : totalCount;
  ScratchSpace<T> tmpScratch(tmpScratchSize, FUN_, "tmpStorage");

  // Pack strided data
  const T* data = 0;
  if(memoryStride == 1) {
    data = array;
  } else {
    for(MPI_Offset i = 0; i < totalCount; i++) {
      tmpScratch[i] = array[memoryStride * i];
    }
    data = tmpScratch.data();
  }

  // labels:IO this is a bugfix for the case when the last process in the communicator
  // has zero elements to write, which resulted in the error:
  //"NetCDF returns status -40: Index exceeds dimension bound"
  // solution taken from here:
  // http://lists.mcs.anl.gov/pipermail/parallel-netcdf/2004-December/000388.html
  if(count[0] == 0) {
    start[0] = 0;
  }

  // Write array
  if(noChunks == 1) {
    // If number of chunks is one, write everything at once
    if(diskStride == 1) {
      status = pnetcdf_traits<T>::ncmpi_put_vara_type_all(b_ncId, varId, start, count, data);
    } else {
      status = pnetcdf_traits<T>::ncmpi_put_vars_type_all(b_ncId, varId, start, count, &diskStride, data);
    }
    b_error(status, name, AT_);
  } else {
    // Write in chunks
    MPI_Offset chunkSize = count[0] / noChunks;
    if(count[0] % noChunks > 0) {
      chunkSize += 1;
    }

    // Determine number of entries for a fixed first dimension index
    size_type nDSize = 1;
    for(size_type d = 1; d < noDims; d++) {
      nDSize *= count[d];
    }

    ScratchSpace<MPI_Offset> start_(noDims, FUN_, "start_");
    ScratchSpace<MPI_Offset> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = max(min(chunkSize, count[0] - i * chunkSize), 0ll);
      const T* data_ = data + min(i * chunkSize * nDSize, (count[0] - 1) * nDSize);
      if(diskStride == 1) {
        status = pnetcdf_traits<T>::ncmpi_put_vara_type_all(b_ncId, varId, &start_[0], &count_[0], data_);
      } else {
        status = pnetcdf_traits<T>::ncmpi_put_vars_type_all(b_ncId, varId, &start_[0], &count_[0], &diskStride, data_);
      }
      b_error(status, name, AT_);
    }
  }
}


/**
 * \brief Writes array data to file (string version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] array Pointer to data in memory.
 * \param[in] name Name of the dataset.
 * \param[in] noDims Number of array dimensions
 * \param[in] start Start offset in the *file* (i.e. the domain should start
 *                  writing from here).
 * \param[in] count Number of elements written from this domain.
 * \param[in] memoryStride Stride in *memory* between data values.
 */
template <>
void ParallelIoPNetcdf::b_writeArray(const MString* array, const MString& name, const size_type noDims,
                                     MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                     const size_type noChunks, MPI_Offset diskStride) {
  TRACE();

  b_ncEndDef();

  // Get variable id
  MInt varId;
  MInt status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine total number of strings (last dimension is string length)
  size_type totalCount = 1;
  for(size_type d = 0; d < noDims - 1; d++) {
    totalCount *= count[d];
  }

  // Determine length of one string
  size_type strLen = count[noDims - 1];

  // Create temporary storage space if needed and set data pointers
  MInt tmpScratchSize = (memoryStride == 1 ? 1 : totalCount);
  ScratchSpace<MString> tmpScratch(tmpScratchSize, FUN_, "tmpStorage");

  // Pack strided data
  const MString* data = 0;
  if(memoryStride == 1) {
    data = array;
  } else {
    for(MPI_Offset i = 0; i < totalCount; i++) {
      tmpScratch[i] = array[memoryStride * i];
    }
    data = tmpScratch.data();
  }

  // labels:IO this is a bugfix for the case when the last process in the communicator
  // has zero elements to write, which resulted in the error:
  //"NetCDF returns status -40: Index exceeds dimension bound"
  // solution taken from here:
  // http://lists.mcs.anl.gov/pipermail/parallel-netcdf/2004-December/000388.html
  if(count[0] == 0) {
    start[0] = 0;
  }

  // Create buffer for writing
  size_type nCount = totalCount * strLen;
  ScratchSpace<char> buf(nCount, FUN_, "buf");

  // Copy data to buffer
  for(MPI_Offset i = 0; i < totalCount; i++) {
    data[i].copy(&buf[i * strLen], strLen, 0);
  }

  // Write array
  if(noChunks == 1) {
    // If number of chunks is one, write everything at once
    if(diskStride == 1) {
      status = ncmpi_put_vara_text_all(b_ncId, varId, start, count, &buf[0]);
    } else {
      status = ncmpi_put_vars_text_all(b_ncId, varId, start, count, &diskStride, &buf[0]);
    }
    b_error(status, name, AT_);
  } else {
    // Write in chunks
    MPI_Offset chunkSize = count[0] / noChunks;
    if(count[0] % noChunks > 0) {
      chunkSize += 1;
    }

    // Determine number of entries for a fixed first dimension index
    size_type nDSize = 1;
    for(size_type d = 1; d < noDims; d++) {
      nDSize *= count[d];
    }

    ScratchSpace<MPI_Offset> start_(noDims, FUN_, "start_");
    ScratchSpace<MPI_Offset> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = max(min(chunkSize, count[0] - i * chunkSize), 0ll);
      const char* buf_ = &buf[0] + min(i * chunkSize * nDSize, (count[0] - 1) * nDSize);
      if(diskStride == 1) {
        status = ncmpi_put_vara_text_all(b_ncId, varId, &start_[0], &count_[0], &buf_[0]);
      } else {
        status = ncmpi_put_vars_text_all(b_ncId, varId, &start_[0], &count_[0], &diskStride, &buf_[0]);
      }
      b_error(status, name, AT_);
    }
  }
}


/**
 * \brief Writes scalar data to file (generic version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] scalar Value to write.
 * \param[in] name Name of the dataset.
 */
template <class T>
void ParallelIoPNetcdf::b_writeScalar(const T scalar, const MString& name) {
  TRACE();

  b_ncEndDef();
  MInt status;

  // Get variable id
  MInt varId;
  status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine offsets
  MPI_Offset start = 0;
  MPI_Offset count = 1;

  // Write scalar
  status = pnetcdf_traits<T>::ncmpi_put_vara_type_all(b_ncId, varId, &start, &count, &scalar);
  b_error(status, name, AT_);
}


/**
 * \brief Writes scalar data to file (string version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] scalar Value to write.
 * \param[in] name Name of the dataset.
 */
template <>
void ParallelIoPNetcdf::b_writeScalar(const MString& scalar, const MString& name) {
  TRACE();

  b_ncEndDef();
  MInt status;

  // Get variable id
  MInt varId;
  status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine offsets
  MPI_Offset start = 0;
  MPI_Offset count = 1;

  // Write scalar
  status = ncmpi_put_vara_text_all(b_ncId, varId, &start, &count, scalar.c_str());
  b_error(status, name, AT_);
}

template <class T>
void ParallelIoPNetcdf::b_readArray(T* array, const MString& path, const MString& name, const size_type noDims,
                                    const size_type* start, const size_type* count) {
  TRACE();

  (void)array;
  (void)path;
  (void)name;
  (void)noDims;
  (void)start;
  (void)count;
  mTerm(1, AT_, "Group functionality not supported by PNetcdf backend");
}


/**
 * \brief Read array data from file (generic version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] array Pointer to data in memory.
 * \param[in] name Name of the dataset.
 * \param[in] noDims Number of array dimensions
 * \param[in] start Start offset in the *file* (i.e. the domain should start
 *                  reading from here).
 * \param[in] count Number of elements read from this domain.
 * \param[in] memoryStride Stride in *memory* between data values.
 *
 * \details This method is the exact counterpart to h5WriteArray().
 */
template <class T>
void ParallelIoPNetcdf::b_readArray(T* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                    MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                    MPI_Offset diskStride) {
  TRACE();

  b_ncEndDef();

  // Get variable id
  MInt varId;
  MInt status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine total data count
  size_type totalCount = 1;
  for(size_type d = 0; d < noDims; d++) {
    totalCount *= count[d];
  }

  // Create temporary storage space if needed and set data pointers
  MInt tmpScratchSize = (memoryStride == 1 ? 1 : totalCount);
  ScratchSpace<T> tmpScratch(tmpScratchSize, FUN_, "tmpStorage");
  T* data = nullptr;
  if(memoryStride == 1) {
    data = array;
  } else {
    data = tmpScratch.data();
  }

  // See b_writeArray(T*) for explanation
  if(count[0] == 0) {
    start[0] = 0;
  }

  // Read array
  if(noChunks == 1) {
    // If number of chunks is one, read everything at once
    if(diskStride == 1) {
      status = pnetcdf_traits<T>::ncmpi_get_vara_type_all(b_ncId, varId, start, count, data);
    } else {
      status = pnetcdf_traits<T>::ncmpi_get_vars_type_all(b_ncId, varId, start, count, &diskStride, data);
    }
    b_error(status, name, AT_);
  } else {
    // Read in chunks
    MPI_Offset chunkSize = count[0] / noChunks;
    if(count[0] % noChunks > 0) {
      chunkSize += 1;
    }

    // Determine number of entries for a fixed first dimension index
    size_type nDSize = 1;
    for(size_type d = 1; d < noDims; d++) {
      nDSize *= count[d];
    }

    ScratchSpace<MPI_Offset> start_(noDims, FUN_, "start_");
    ScratchSpace<MPI_Offset> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = max(min(chunkSize, count[0] - i * chunkSize), 0ll);
      T* data_ = data + min(i * chunkSize * nDSize, (count[0] - 1) * nDSize);
      if(diskStride == 1) {
        status = pnetcdf_traits<T>::ncmpi_get_vara_type_all(b_ncId, varId, &start_[0], &count_[0], data_);
      } else {
        status = pnetcdf_traits<T>::ncmpi_get_vars_type_all(b_ncId, varId, &start_[0], &count_[0], &diskStride, data_);
      }
      b_error(status, name, AT_);
    }
  }

  // Unpack strided data if necessary
  if(memoryStride != 1) {
    for(MPI_Offset i = 0; i < totalCount; i++) {
      array[memoryStride * i] = tmpScratch[i];
    }
  }
}


/**
 * \brief Read array data from file (string version). <b>[MPI]</b>
 *
 * \author Ansgar Niemoeller, Konstantin Froehlich
 * \date 2015-08-19
 *
 * \param[in] array Pointer to data in memory.
 * \param[in] name Name of the dataset.
 * \param[in] noDims Number of array dimensions
 * \param[in] start Start offset in the *file* (i.e. the domain should start
 *                  reading from here).
 * \param[in] count Number of elements read from this domain.
 * \param[in] memoryStride Stride in *memory* between data values.
 *
 * \details This method is the exact counterpart to h5WriteArray().
 */
template <>
void ParallelIoPNetcdf::b_readArray(MString* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                    MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                    MPI_Offset diskStride) {
  TRACE();

  b_ncEndDef();

  // Get variable id
  MInt varId;
  MInt status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine total data count
  size_type totalCount = 1;
  for(size_type d = 0; d < noDims - 1; d++) {
    totalCount *= count[d];
  }

  // Determine length of one string
  size_type strLen = count[noDims - 1];

  // Create temporary storage space if needed and set data pointers
  MInt tmpScratchSize = (memoryStride == 1 ? 1 : totalCount);
  ScratchSpace<MString> tmpScratch(tmpScratchSize, FUN_, "tmpStorage");
  MString* data = nullptr;
  if(memoryStride == 1) {
    data = array;
  } else {
    data = tmpScratch.data();
  }

  // Create buffer for reading
  size_type nCount = totalCount * strLen;
  ScratchSpace<char> buf(nCount, FUN_, "buf");

  // Read array
  if(noChunks == 1) {
    // If number of chunks is one, read everything at once
    if(diskStride == 1) {
      status = ncmpi_get_vara_text_all(b_ncId, varId, start, count, &buf[0]);
    } else {
      status = ncmpi_get_vars_text_all(b_ncId, varId, start, count, &diskStride, &buf[0]);
    }
    b_error(status, name, AT_);

    // Extract strings from buffer
    for(size_type i = 0; i < totalCount; i++) {
      MString tmp;
      tmp.append(&buf[i * strLen], strLen);
      data[i].append(tmp.c_str(), 0, strLen);
    }
  } else {
    // Read in chunks
    MPI_Offset chunkSize = count[0] / noChunks;
    if(count[0] % noChunks > 0) {
      chunkSize += 1;
    }

    ScratchSpace<MPI_Offset> start_(noDims, FUN_, "start_");
    ScratchSpace<MPI_Offset> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = max(min(chunkSize, count[0] - i * chunkSize), 0ll);
      if(diskStride == 1) {
        status = ncmpi_get_vara_text_all(b_ncId, varId, &start_[0], &count_[0], &buf[0]);
      } else {
        status = ncmpi_get_vars_text_all(b_ncId, varId, &start_[0], &count_[0], &diskStride, &buf[0]);
      }

      b_error(status, name, AT_);

      // Extract strings from buffer
      for(size_type j = start_[0]; j < totalCount; j++) {
        MString tmp;
        tmp.append(&buf[(j - start_[0]) * strLen], strLen);
        data[j].append(tmp, 0, strLen);
      }
    }
  }

  // Unpack strided data if necessary
  if(memoryStride != 1) {
    for(MPI_Offset i = 0; i < totalCount; i++) {
      array[memoryStride * i] = tmpScratch[i];
    }
  }
}


/**
 * \brief Read scalar data from file (generic version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[out] scalar Value to read.
 * \param[in] name Name of the dataset.
 *
 * \details This method is the exact counterpart to h5WriteScalar().
 */
template <class T>
void ParallelIoPNetcdf::b_readScalar(T* scalar, const MString& name) {
  TRACE();

  b_ncEndDef();
  MInt status;

  // Get variable id
  MInt varId;
  status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine offsets
  MPI_Offset start = 0;
  MPI_Offset count = 1;

  // Read scalar
  status = pnetcdf_traits<T>::ncmpi_get_vara_type_all(b_ncId, varId, &start, &count, scalar);
  b_error(status, name, AT_);
}


/**
 * \brief Read scalar data from file (float version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[out] scalar Value to read.
 * \param[in] name Name of the dataset.
 *
 * \details This method is the exact counterpart to h5WriteScalar().
 */
template <>
void ParallelIoPNetcdf::b_readScalar(MString* scalar, const MString& name) {
  TRACE();

  b_ncEndDef();
  MInt status;

  // Get variable id
  MInt varId;
  status = ncmpi_inq_varid(b_ncId, name.c_str(), &varId);
  b_error(status, name, AT_);

  // Determine offsets
  MPI_Offset start = 0;
  MPI_Offset count = 1;

  // Read scalar (one char)
  status = ncmpi_get_vara_text_all(b_ncId, varId, &start, &count, (char*)scalar);
  b_error(status, name, AT_);
}

//------------------------------------------------------------------------------
// Attribute methods
//------------------------------------------------------------------------------


/**
 * \brief Set an attribute in the file (generic version).
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar1
 *
 * \param[in] value Attribute value.
 * \param[in] name Attribute name.
 * \param[in] datasetName Dataset (scalar or array) to which the attribute is
 *                        attached (leave empty for file attributes).
 */
template <class T>
void ParallelIoPNetcdf::b_setAttribute(const T* value, const MString& name, const MString& datasetName,
                                       const size_type totalCount) {
  TRACE();

  // Determine variable id
  MInt varId;
  if(datasetName.empty()) {
    varId = NC_GLOBAL;
  } else {
    MInt status = ncmpi_inq_varid(b_ncId, datasetName.c_str(), &varId);
    b_error(status, datasetName, AT_);
  }

  // If attribute does not exist, go to define mode
  if(!b_hasAttribute(name, datasetName)) {
    b_ncRedef();
  }

  // Write attribute
  MInt status =
      pnetcdf_traits<T>::ncmpi_put_att_type(b_ncId, varId, name.c_str(), pnetcdf_traits<T>::type(), totalCount, value);
  b_error(status, datasetName + "::" + name, AT_);
}


/**
 * \brief Set an attribute in the file (string version).
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar
 *
 * \param[in] value Attribute value.
 * \param[in] name Attribute name.
 * \param[in] datasetName Dataset (scalar or array) to which the attribute is
 *                        attached (leave empty for file attributes).
 */
template <>
void ParallelIoPNetcdf::b_setAttribute(const MString* value, const MString& name, const MString& datasetName,
                                       const size_type totalCount) {
  TRACE();

  if(totalCount > 1) {
    mTerm(1, AT_, "Array of strings attributes not yet supported.");
  }

  // Determine variable id
  MInt varId;
  if(datasetName.empty()) {
    varId = NC_GLOBAL;
  } else {
    // If this is a data file, get data file-specific name
    MInt status = ncmpi_inq_varid(b_ncId, datasetName.c_str(), &varId);
    b_error(status, datasetName, AT_);
  }

  // If attribute does not exist or is of greater size, go to define mode
  if(!b_hasAttribute(name, datasetName)) {
    b_ncRedef();
  } else {
    MPI_Offset length;
    MInt status = ncmpi_inq_attlen(b_ncId, varId, name.c_str(), &length);
    b_error(status, datasetName + "::" + name, AT_);

    if(length < static_cast<MPI_Offset>(value->length())) {
      b_ncRedef();
    }
  }

  // Write attribute
  MInt status = ncmpi_put_att_text(b_ncId, varId, name.c_str(), value->length(), (*value).c_str());
  b_error(status, datasetName + "::" + name, AT_);
}


/**
 * \brief Retrieve an attribute from file (generic version).
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar
 *
 * \param[out] value Attribute value
 * \param[in] name Attribute name
 * \param[in] datasetName Dataset (scalar or array) to which the attribute is
 *                        attached (leave empty for file attributes).
 */
template <class T>
void ParallelIoPNetcdf::b_getAttribute(T* const value, const MString& name, const MString& datasetName,
                                       const size_type totalCount) {
  TRACE();

  // Determine variable id
  MInt varId;
  if(datasetName.empty()) {
    varId = NC_GLOBAL;
  } else {
    MInt status = ncmpi_inq_varid(b_ncId, datasetName.c_str(), &varId);
    b_error(status, datasetName, AT_);
  }

  // Get attribute length
  MInt status;
  MPI_Offset length;
  status = ncmpi_inq_attlen(b_ncId, varId, name.c_str(), &length);
  b_error(status, datasetName + "::" + name, AT_);

  if(length != (MPI_Offset)totalCount) {
    TERMM(1, "Requested attribute (" + name + ") has different number of elements (" + to_string(length)
                 + ") than specified (" + to_string(totalCount)
                 + "). Use getAttributeCount() to query number of elements first");
  }

  // Read attribute
  status = pnetcdf_traits<T>::ncmpi_get_att_type(b_ncId, varId, name.c_str(), value);
  b_error(status, datasetName + "::" + name, AT_);
}


/**
 * \brief Retrieve an attribute from file (string version).
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar
 *
 * \param[out] value Attribute value
 * \param[in] name Attribute name
 * \param[in] datasetName Dataset (scalar or array) to which the attribute is
 *                        attached (leave empty for file attributes).
 */
template <>
void ParallelIoPNetcdf::b_getAttribute(MString* const value, const MString& name, const MString& datasetName,
                                       const size_type totalCount) {
  TRACE();

  if(totalCount > 1) {
    mTerm(1, AT_, "Array of strings attributes not yet supported.");
  }

  // Determine variable id
  MInt varId;
  if(datasetName.empty()) {
    varId = NC_GLOBAL;
  } else {
    MInt status = ncmpi_inq_varid(b_ncId, datasetName.c_str(), &varId);
    b_error(status, datasetName, AT_);
  }

  // Get attribute length
  MInt status;
  MPI_Offset length;
  status = ncmpi_inq_attlen(b_ncId, varId, name.c_str(), &length);
  b_error(status, datasetName + "::" + name, AT_);

  // Read attribute
  ScratchSpace<MChar> tmpScratch(length, FUN_, "tmpScratch");
  status = ncmpi_get_att_text(b_ncId, varId, name.c_str(), tmpScratch.data());
  b_error(status, datasetName + "::" + name, AT_);
  value->assign(tmpScratch.data(), length);
}


void ParallelIoPNetcdf::b_getAttributeCount(const MString& name, size_type* totalCount, const MString& datasetName) {
  TRACE();

  // Determine variable id
  MInt varId;
  if(datasetName.empty()) {
    varId = NC_GLOBAL;
  } else {
    MInt status = ncmpi_inq_varid(b_ncId, datasetName.c_str(), &varId);
    b_error(status, datasetName, AT_);
  }

  // Get attribute length
  MInt status;
  MPI_Offset length;
  status = ncmpi_inq_attlen(b_ncId, varId, name.c_str(), &length);
  b_error(status, datasetName + "::" + name, AT_);

  *totalCount = (size_type)length;
}


//------------------------------------------------------------------------------
// Auxiliary methods
//------------------------------------------------------------------------------
/**
 * \brief Check the status code of a HDF5 operation and output a meaningful
 *        message.
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar
 *
 * \param[in] status      HDF5 status code (negative if unexspected behaviour)
 * \param[in] name        Name of the file / variable / attribute
 * \param[in] location    Function which called the error message
 *
 * \return The status message if there was an error, or zero if everything is
 *  ok.
 */
void ParallelIoPNetcdf::b_error(MInt status, const MString& name, const MString& location) {
  // TRACE(); //<- this function is called to often and not really important.
  if(status != NC_NOERR) {
    cerr << endl;
    cerr << "*** ERROR in parallelio_pnetcdf ***" << endl;
    cerr << "NetCDF error in '" << location << "'" << endl;
    cerr << "NetCDF returns status " << status << ": " << ncmpi_strerror(status) << endl;
    cerr << "The file/variable/attribute in question was: " << name << endl;
    cerr << endl;
    TERMM(1, "NetCDF error in ParallelIo.");
  }
}


// Explicit instantiations for all supported types
// Write methods
template void ParallelIoPNetcdf::b_writeArray(const MFloat* const array, const MString& name, const size_type noDims,
                                              MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                              const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_writeArray(const MInt* const array, const MString& name, const size_type noDims,
                                              MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                              const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_writeArray(const MLong* const array, const MString& name, const size_type noDims,
                                              MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                              const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_writeArray(const MUchar* const array, const MString& name, const size_type noDims,
                                              MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                              const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_writeScalar(const MFloat scalar, const MString& name);
template void ParallelIoPNetcdf::b_writeScalar(const MInt scalar, const MString& name);
template void ParallelIoPNetcdf::b_writeScalar(const MLong scalar, const MString& name);
template void ParallelIoPNetcdf::b_writeScalar(const MUchar scalar, const MString& name);
// Read methods
template void ParallelIoPNetcdf::b_readArray(MFloat* const array, const MString& name, const size_type noDims,
                                             MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                             const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_readArray(MInt* const array, const MString& name, const size_type noDims,
                                             MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                             const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_readArray(MLong* const array, const MString& name, const size_type noDims,
                                             MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                             const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_readArray(MUchar* const array, const MString& name, const size_type noDims,
                                             MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                             const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoPNetcdf::b_readScalar(MFloat* scalar, const MString& name);
template void ParallelIoPNetcdf::b_readScalar(MInt* scalar, const MString& name);
template void ParallelIoPNetcdf::b_readScalar(MLong* scalar, const MString& name);
template void ParallelIoPNetcdf::b_readScalar(MUchar* scalar, const MString& name);
// Attribute methods
template void ParallelIoPNetcdf::b_setAttribute(const MFloat* value, const MString& name, const MString& datasetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_setAttribute(const MInt* value, const MString& name, const MString& datasetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_setAttribute(const MLong* value, const MString& name, const MString& datasetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_setAttribute(const MUchar* value, const MString& name, const MString& datasetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_setAttribute(const MUlong* value, const MString& name, const MString& datasetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_getAttribute(MFloat* const value, const MString& name, const MString& datagetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_getAttribute(MInt* const value, const MString& name, const MString& datagetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_getAttribute(MLong* const value, const MString& name, const MString& datagetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_getAttribute(MUchar* const value, const MString& name, const MString& datagetName,
                                                const size_type totalCount);
template void ParallelIoPNetcdf::b_getAttribute(MUlong* const value, const MString& name, const MString& datagetName,
                                                const size_type totalCount);
