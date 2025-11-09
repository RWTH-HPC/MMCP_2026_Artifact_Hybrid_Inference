// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "parallelio.h"

#if defined(WITH_HDF5)
#include "parallelio_hdf5.h"

// Due to consistency with Netcdf: Maximal name length for datasets is set
// to NC_MAX_NAME
#include "hdf5_hl.h"


#include <cstdlib>
#include <cstring>
#include <numeric>
#include <string.h>
#include <sys/stat.h>
#include "COMM/mpioverride.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "typetraits.h"

#if not defined(MAIA_WINDOWS)
#include <pwd.h>
#include <unistd.h>
#include "maiapnetcdf.h"
#else
#define NC_MAX_NAME 256
#include <WinSock2.h>
#include <direct.h>
#endif

using namespace maia;
using namespace parallel_io;
using namespace std;

//------------------------------------------------------------------------------
// HDF5-specific methods
//------------------------------------------------------------------------------

// Use unnamed namespace for some HDF5-specific magic (a.k.a. type traits)
namespace {

template <class DataType>
struct hdf5_traits {};

// MFloat traits
template <>
struct hdf5_traits<MFloat> {
  // Corresponding HDF5 data type
  static hid_t type() { return H5T_NATIVE_DOUBLE; }
};

// MInt traits
template <>
struct hdf5_traits<MInt> {
  // Corresponding HDF5 data type
  static hid_t type() { return H5T_NATIVE_INT; }
};

// MLong traits
template <>
struct hdf5_traits<MLong> {
  // Corresponding HDF5 data type
  static hid_t type() { return H5T_NATIVE_LONG; }
};

// MChar traits
template <>
struct hdf5_traits<MChar> {
  // Corresponding HDF5 data type
  static hid_t type() { return H5T_NATIVE_CHAR; }
};

// MUchar traits
template <>
struct hdf5_traits<MUchar> {
  // Corresponding HDF5 data type
  static hid_t type() { return H5T_NATIVE_UCHAR; }
};

// MUlong traits
template <>
struct hdf5_traits<MUlong> {
  // Corresponding HDF5 data type
  static hid_t type() { return H5T_NATIVE_ULLONG; }
};

} // namespace

/*
// HDF5 specific IO is not tested with valgrind, yet. "maiapnetcdf.h" handles
// some issues regarding valgrind and errors of the third party library Parallel
// Netcdf (-> wiki (Version: Nov 2015)).
// There is no such "hdf5.h" for Hdf5. Additionally some testcases might
// suppress valgrind errors if you run them with Parallel Netcdf. Probably,
// these testcases have also surpress according Hdf5 errors.
//
// The HDF5 API defines aliases for datatypes which are used for specific
// purposes. Examples are:
// hid_t: type for "managing references to nodes"(identifier). Each node
// reference is represented as an integer. All identifiers have to be closed
// after usage!
// herr_t: type for handling error codes. Mainly used to signal whether a
// function call was successful.
// hsize_t: represents a native multiple-precision integer.
*/

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
 * \param[in] Not used by Hdf5 specific version.
 *
 */
MBool ParallelIoHdf5::b_isValidFile(const MString& name, const MPI_Comm& /*not used*/) {
  TRACE();

  // WARNING: This method is untested, yet. It is implemented due to be
  // consistent with parallel NetCDF. However, it was not used by any of the
  // testcases.If your code uses this part of the code, please make sure that
  // the I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().
  MBool returnValue;

  herr_t status;
  hid_t fileId;
  fileId = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(fileId > 0) {
    returnValue = true;
    status = H5Fclose(fileId);
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
const MString ParallelIoHdf5::b_fileExt() {
  TRACE();

  return ".Hdf5";
}


//------------------------------------------------------------------------------
// Constructor & Destructor
//------------------------------------------------------------------------------
/**
 * \brief Creates a new object to read and write *big* data files. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-2
 * \author (Modified by: )Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-03-01
 *
 * \param[in] fileName
 * \param[in] fileMode The file mode can be either maia::parallel_io::PIO_CREATE,
 *                     maia::parallel_io::PIO_APPEND, or maia::parallel_io::PIO_READ.
 * \param[in] mpiComm The MPI communicator that should be used to open/create
 *                    the file.
 */
ParallelIoHdf5::ParallelIoHdf5(const MString& fileName, MInt fileMode, const MPI_Comm& mpiComm)
  : ParallelIoBase<ParallelIoHdf5>(fileName, fileMode, mpiComm) {
  TRACE();

#ifdef DISABLE_OUTPUT
  if(m_fileMode != PIO_READ) return;
#endif

  // Store MPI IO communicator
  b_h5FileXferHandle = H5Pcreate(H5P_FILE_ACCESS);
  herr_t status = H5Pset_fapl_mpio(b_h5FileXferHandle, m_mpiComm, m_mpiInfo);
  b_error(status, m_fileName, AT_);

  // Setting b_h5DatasetXferHandle to collective access
  b_h5DatasetXferHandle = H5Pcreate(H5P_DATASET_XFER);

  status = H5Pset_dxpl_mpio(b_h5DatasetXferHandle, H5FD_MPIO_COLLECTIVE);
  b_error(status, m_fileName, AT_);

  switch(m_fileMode) {
    case PIO_CREATE: {
      b_h5Id = H5Fcreate(m_fileName.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, b_h5FileXferHandle);
      b_error(b_h5Id, m_fileName + " :: CREATE", AT_);
    } break;

    case PIO_REPLACE: {
      b_h5Id = H5Fcreate(m_fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, b_h5FileXferHandle);
      b_error(b_h5Id, m_fileName + " :: REPLACE", AT_);
    } break;

    case PIO_APPEND: {
      // Attempt to open an existing file to append data
      b_h5Id = H5Fopen(m_fileName.c_str(), H5F_ACC_RDWR, b_h5FileXferHandle);
      b_error(b_h5Id, m_fileName + " :: APPEND", AT_);
    } break;

    case PIO_READ: {
      b_h5Id = H5Fopen(m_fileName.c_str(), H5F_ACC_RDONLY, b_h5FileXferHandle);
      b_error(b_h5Id, m_fileName + " :: READ", AT_);
    } break;

    default: {
      mTerm(1, AT_, "Unsupported file mode.");
    } break;
  }
}

/**
 * \brief Close open identifiers and release memory. <b>[MPI]</b>
 *
 * \author Konstantin Froehlich
 * \date 2015-07-01
 */
ParallelIoHdf5::~ParallelIoHdf5() {
  TRACE();

  // Close all identifiers used for collective writing
  herr_t status = H5Pclose(b_h5DatasetXferHandle);
  b_error(status, m_fileName, AT_);
  status = H5Pclose(b_h5FileXferHandle);
  b_error(status, m_fileName, AT_);
  status = H5Fclose(b_h5Id);
  b_error(status, m_fileName, AT_);

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
 * \brief Adds all additional header information that are needed in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \details This method must be called exactly once for newly created files.
 */
void ParallelIoHdf5::b_saveHeader() {
  TRACE();

  b_addAdditionalHeader();
  b_writeAdditionalData();
}


// Forward declaration of specialization for use in b_addAdditionalHeader
template <>
void ParallelIoHdf5::b_setAttribute(const MString* value, const MString& name, const MString& datasetName,
                                    const size_type totalCount);


/**
 * \brief Write additional headers to file (e.g. grid file name, creation date
 *        etc.). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-12
 */
void ParallelIoHdf5::b_addAdditionalHeader() {
  TRACE();

  // For a newly created or modified file, add some meta information
  if(m_fileMode != PIO_READ) {
    const MInt maxNoChars = 256;

    // Get all meta-data on root process & communicate, since some of the C
    // functions might not be thread-safe (at least getpwuid() is not)
    MChar user[maxNoChars];
    MChar host[maxNoChars];
    MChar dir[maxNoChars];
    MChar date[maxNoChars];

    // Create object & initialize data
    fill(user, user + maxNoChars, '\0');
    fill(host, host + maxNoChars, '\0');
    fill(dir, dir + maxNoChars, '\0');
    fill(date, date + maxNoChars, '\0');

    if(m_domainId == 0) {
      // Gets the current username
      passwd* p;
#if not defined(MAIA_WINDOWS)
      p = getpwuid(getuid());
      if(p) {
        strncpy(user, p->pw_name, maxNoChars - 1);
      } else {
        strncpy(user, "n/a", maxNoChars - 1);
      }
#else
      strncpy(user, "windows", maxNoChars - 1);
#endif

      // Gets the current hostname
      gethostname(host, maxNoChars - 1);

      // Gets the current directory
#if defined(MAIA_WINDOWS)
      _getcwd(dir, maxNoChars - 1);
#else
      char* retv = getcwd(dir, maxNoChars - 1);
      ASSERT(retv != nullptr, "");
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
    const MInt noItems = 4;
    MChar buffer[noItems * maxNoChars];
    memcpy(buffer + 0 * maxNoChars, user, maxNoChars);
    memcpy(buffer + 1 * maxNoChars, host, maxNoChars);
    memcpy(buffer + 2 * maxNoChars, dir, maxNoChars);
    memcpy(buffer + 3 * maxNoChars, date, maxNoChars);

    // Broadcast time from rank 0 to ensure that every rank has the same
    // information
    MPI_Bcast(&buffer, 4 * maxNoChars, MPI_CHAR, 0, m_mpiComm, AT_, "buffer");

    // Unpack data
    memcpy(user, buffer + 0 * maxNoChars, maxNoChars);
    memcpy(host, buffer + 1 * maxNoChars, maxNoChars);
    memcpy(dir, buffer + 2 * maxNoChars, maxNoChars);
    memcpy(date, buffer + 3 * maxNoChars, maxNoChars);

    MString version = MString(XSTRINGIFY(MAIA_VERSION_STRING));
    MString build = MString(XSTRINGIFY(MAIA_COMPILER_STRING)) + " " + MString(XSTRINGIFY(MAIA_BUILD_TYPE_STRING)) + " ("
                    + MString(XSTRINGIFY(MAIA_COMPILER_VERSION_STRING)) + ")";

    if(m_fileMode == PIO_CREATE || m_fileMode == PIO_REPLACE) {
      // Add file attributes only needed for creation
      ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(user), "meta_creation_user");
      ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(host), "meta_creation_host");
      ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(dir), "meta_creation_directory");
      ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(date), "meta_creation_date");
      ParallelIoBase<ParallelIoHdf5>::setAttribute(m_noDomains, "meta_creation_noDomains");
      ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(date), "_meta_creation_date");
      ParallelIoBase<ParallelIoHdf5>::setAttribute(version, "_meta_creation_revision");
      ParallelIoBase<ParallelIoHdf5>::setAttribute(build, "_meta_creation_build");
    }

    // Add file attributes that should be set at each modification
    ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(user), "meta_lastModified_user");
    ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(host), "meta_lastModified_host");
    ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(dir), "meta_lastModified_directory");
    ParallelIoBase<ParallelIoHdf5>::setAttribute(MString(date), "meta_lastModified_date");
    ParallelIoBase<ParallelIoHdf5>::setAttribute(m_noDomains, "meta_lastModified_noDomains");
    ParallelIoBase<ParallelIoHdf5>::setAttribute(version, "_meta_lastModified_revision");
    ParallelIoBase<ParallelIoHdf5>::setAttribute(build, "_meta_lastModified_build");
  }
}


/**
 * \brief Write additional data to file. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-26
 */
void ParallelIoHdf5::b_writeAdditionalData() {
  TRACE();

  // At the moment, nothing happens here
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
 *                 maia::parallel_io::PIO_STRING maia::parallel_io::PIO_UCHAR).
 * \param[in] name Name of the array (must not be empty).
 * \param[in] noDims Number of array dimensions
 * \param[in] totalCount Total size of the array in each dimension.
 *
 */
void ParallelIoHdf5::b_defineArray(maiabd_type type, const MString& name, size_type noDims, size_type* totalCount) {
  TRACE();

  // Choose the HDF5 specific native datatype or a string of variable length
  herr_t status = -1;
  hid_t dtype_id;
  switch(type) {
    case PIO_FLOAT: {
      dtype_id = H5T_NATIVE_DOUBLE;
      break;
    }
    case PIO_INT: {
      dtype_id = H5T_NATIVE_INT;
      break;
    }
    case PIO_LONG: {
      dtype_id = H5T_NATIVE_LONG;
      break;
    }
    case PIO_STRING: {
      dtype_id = H5Tcopy(H5T_C_S1);
      b_error(dtype_id, name, AT_);
      status = H5Tset_size(dtype_id, H5T_VARIABLE);
      b_error(dtype_id, name, AT_);
      break;
    }
    case PIO_UCHAR: {
      dtype_id = H5T_NATIVE_UCHAR;
      break;
    }
    case PIO_ULONGLONG: {
      dtype_id = H5T_NATIVE_ULLONG;
      break;
    }

    default: {
      TERMM(1, "Invalid ParallelIo data type!");
    }
  }

  // Allocate Dataspace as big as specified by totalCount and noDims
  hid_t dspace_id = H5Screate_simple(noDims, (hsize_t*)(totalCount), nullptr);
  b_error(dspace_id, name, AT_);
  hid_t dset_id = H5Dcreate2(b_h5Id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  b_error(dset_id, name, AT_);

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  if(type == PIO_STRING) {
    status = H5Tclose(dtype_id);
    b_error(status, name, AT_);
  }
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
void ParallelIoHdf5::b_defineScalar(maiabd_type type, const MString& name) {
  TRACE();

  // Choose the HDF5 specific native datatype or a string of variable length
  herr_t status = -1;
  hid_t dtype_id;
  switch(type) {
    case PIO_FLOAT: {
      dtype_id = H5T_NATIVE_DOUBLE;
      break;
    }
    case PIO_INT: {
      dtype_id = H5T_NATIVE_INT;
      break;
    }
    case PIO_LONG: {
      dtype_id = H5T_NATIVE_LONG;
      break;
    }
    case PIO_STRING: {
      dtype_id = H5Tcopy(H5T_C_S1);
      b_error(dtype_id, name, AT_);
      status = H5Tset_size(dtype_id, H5T_VARIABLE);
      b_error(dtype_id, name, AT_);
      break;
    }
    case PIO_UCHAR: {
      dtype_id = H5T_NATIVE_UCHAR;
      break;
    }
    case PIO_ULONGLONG: {
      dtype_id = H5T_NATIVE_ULLONG;
      break;
    }
    default: {
      TERMM(1, "Invalid ParallelIo data type!");
    }
  }

  // Create a dataspace and dataset for a scalar
  hid_t dspace_id = H5Screate(H5S_SCALAR);
  b_error(dspace_id, name, AT_);
  hid_t dset_id = H5Dcreate2(b_h5Id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  b_error(dset_id, name, AT_);

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  if(type == PIO_STRING) {
    status = H5Tclose(dtype_id);
    b_error(status, name, AT_);
  }
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
MBool ParallelIoHdf5::b_hasDataset(const MString& name, const size_type noDimensions) {
  TRACE();

  MBool varExists = false;

  herr_t status = -1;
  hid_t dset_id = -1, dspace_id = -1;
  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);
  htri_t linkExists = H5Lexists(b_h5Id, name.c_str(), link_id);
  b_error(linkExists, name, AT_);
  if(linkExists != 0) {
    dset_id = H5Dopen2(b_h5Id, name.c_str(), H5P_DEFAULT);
    b_error(dset_id, name, AT_);
    dspace_id = H5Dget_space(dset_id);
    b_error(dspace_id, name, AT_);
    int nDims = H5Sget_simple_extent_ndims(dspace_id);
    b_error(nDims, name, AT_);
    if(nDims == noDimensions || noDimensions == -1) {
      varExists = true;
    }
  }

  // Close all previously created identifiers (except b_h5Id)
  if(dset_id > 0) {
    status = H5Sclose(dspace_id);
    b_error(status, name, AT_);
    status = H5Dclose(dset_id);
    b_error(status, name, AT_);
  }
  status = H5Pclose(link_id);
  b_error(status, name, AT_);

  return varExists;
}


/**
 * \brief Check if dataset exists. <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[out]
 * \param[in]
 *
 * \details
 */
MBool ParallelIoHdf5::b_hasDataset(const MString& name, const MString& path) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);
  // Handle whether or not path was given
  MBool exists = false;
  hid_t loc_id = 0;
  MBool checkGroup = false;
  if(path.empty()) {
    loc_id = b_h5Id;
  } else {
    herr_t pathExists = H5Lexists(b_h5Id, path.c_str(), link_id);
    if(pathExists != 0) {
      // Group path exists, opening it
      loc_id = H5Oopen(b_h5Id, path.c_str(), H5P_DEFAULT);
      b_error(loc_id, path, AT_);
    } else {
      // Group path doesn't exist, then object cannot exist
      H5Pclose(link_id);
      return 0;
    }

    checkGroup = true;
  }

  // Check if object exists
  herr_t status = H5Lexists(loc_id, name.c_str(), link_id);
  if(status != 0) {
    exists = true;
  }

  H5Pclose(link_id);
  if(checkGroup) H5Oclose(loc_id);
  return exists;
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
MInt ParallelIoHdf5::b_getDatasetType(const MString& name) {
  TRACE();

  // Determine variable id
  hid_t dset_id = H5Dopen2(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  MInt typeId;
  if(H5Tequal(dtype_id, H5T_NATIVE_INT) > 0) {
    typeId = PIO_INT;
  } else if(H5Tequal(dtype_id, H5T_NATIVE_DOUBLE) > 0) {
    typeId = PIO_FLOAT;
  } else if(H5Tequal(dtype_id, H5T_NATIVE_LONG) > 0) {
    typeId = PIO_LONG;
  } else if(H5Tget_class(dtype_id) == H5T_STRING || H5Tequal(dtype_id, H5T_C_S1) > 0) {
    typeId = PIO_STRING;
  } else if(H5Tequal(dtype_id, H5T_NATIVE_UCHAR) > 0) {
    typeId = PIO_UCHAR;
  } else if(H5Tequal(dtype_id, H5T_NATIVE_ULLONG) > 0) {
    typeId = PIO_ULONGLONG;
  } else {
    typeId = PIO_UNKNOWN_TYPE;
    TERMM(1, "ERROR: Unknown type for dataset " + name + "!");
  }

  // Close all previously created identifiers (except b_h5Id)
  herr_t status = H5Tclose(dtype_id);
  b_error(status, name, AT_);
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);

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
void ParallelIoHdf5::b_getDatasetNames(vector<MString>& names, const size_type dimension) {
  TRACE();

  vector<MString>().swap(names);

  herr_t status = -1;

  // Get access to the root group "/" (currently group functionality of HDF5 is
  // completely ignored, thus it is the only existent group).
  hid_t group_id = H5Gopen(b_h5Id, "/", H5P_DEFAULT);
  b_error(group_id, m_fileName, AT_);
  // Iterate through each object in the group.
  hsize_t noObj = 0;
  status = H5Gget_num_objs(group_id, &noObj);
  b_error(status, m_fileName, AT_);
  for(MInt idx = 0; idx < (signed)noObj; idx++) {
    // Check if the current object is a dataset
    int obj_type = H5Gget_objtype_by_idx(group_id, (size_t)idx);
    b_error(obj_type, m_fileName, AT_);
    if(obj_type == H5G_DATASET) {
      // Determine the name and the dimension of the dataset
      // The maximal name length is NC_MAX_NAME for consistency with PnetCDF
      char DatasetName[NC_MAX_NAME];
      status = H5Gget_objname_by_idx(group_id, (hsize_t)idx, DatasetName, NC_MAX_NAME);
      b_error(status, m_fileName, AT_);
      hid_t dset_id = H5Dopen(b_h5Id, DatasetName, H5P_DEFAULT);
      b_error(dset_id, m_fileName + " :: " + string(DatasetName), AT_);
      hid_t dspace_id = H5Dget_space(dset_id);
      b_error(dspace_id, m_fileName + " :: " + string(DatasetName), AT_);
      int noDims = H5Sget_simple_extent_ndims(dspace_id);
      b_error(noDims, m_fileName + " :: " + string(DatasetName), AT_);
      if(noDims == dimension || dimension == -1) {
        names.emplace_back(DatasetName);
      }
      // Close all previously created identifiers (except b_h5Id)
      status = H5Sclose(dspace_id);
      b_error(dspace_id, m_fileName + " :: " + string(DatasetName), AT_);
      status = H5Dclose(dset_id);
      b_error(dspace_id, m_fileName + " :: " + string(DatasetName), AT_);
    }
  }
  status = H5Gclose(group_id);
  b_error(status, m_fileName, AT_);
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
ParallelIo::size_type ParallelIoHdf5::b_getDatasetNoDims(const MString& name) {
  TRACE();

  // Get variable id
  hid_t dset_id = H5Dopen2(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);

  // Get number of variable dimensions
  hid_t dspace_id = H5Dget_space(dset_id);
  b_error(dspace_id, name, AT_);
  int noDims = H5Sget_simple_extent_ndims(dspace_id);
  b_error(noDims, name, AT_);

  // Close all previously created identifiers (except b_h5Id)
  herr_t status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Dclose(dset_id);
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
ParallelIo::size_type ParallelIoHdf5::b_getDatasetSize(const MString& name, const size_type dimensionId) {
  TRACE();

  // Get variable id
  hid_t dset_id = H5Dopen2(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);

  // Get number of array dimensions
  size_type noDims = b_getDatasetNoDims(name);

  // Get variable size for each dimension
  hid_t dspace_id = H5Dget_space(dset_id);
  b_error(dspace_id, name, AT_);

  ScratchSpace<hsize_t> dimId(noDims, FUN_, "dimId");

  herr_t status = H5Sget_simple_extent_dims(dspace_id, &dimId[0], nullptr);
  b_error(status, name, AT_);
  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  // Close all previously created identifiers (except b_h5Id)
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);

  return static_cast<size_type>(dimId[dimensionId]);
}

//------------------------------------------------------------------------------
// Attribute methods
//------------------------------------------------------------------------------

/**
 * \brief Retrieve an attribute from file (MString version).
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
void ParallelIoHdf5::b_getAttribute(MString* value, const MString& name, const MString& datasetName,
                                    const size_type totalCount) {
  TRACE();

  if(totalCount > 1) {
    mTerm(1, AT_, "String array attributes not supported.");
  }

  // Open the attribute
  hid_t attribute_id;
  if(datasetName == "") { // Attribute is attached to the file
    attribute_id = H5Aopen(b_h5Id, name.c_str(), H5P_DEFAULT);
    b_error(attribute_id, name, AT_);
  } else { // Attribute is attached to a dataset
    attribute_id = H5Aopen_by_name(b_h5Id, datasetName.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    b_error(attribute_id, name, AT_);
  }
  hid_t attribute_type_id = H5Aget_type(attribute_id);
  b_error(attribute_type_id, name, AT_);
  hsize_t length = H5Tget_size(attribute_type_id);

  hid_t filetype = H5Tcopy(H5T_C_S1);
  H5Tset_size(filetype, length);
  ScratchSpace<char> buf(length, FUN_, "buf");

  herr_t status = H5Aread(attribute_id, filetype, &buf[0]);
  b_error(status, name, AT_);

  // Extract strings from buffer
  value[0].append(&buf[0], length);

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Tclose(attribute_type_id);
  b_error(status, name, AT_);
  status = H5Aclose(attribute_id);
  b_error(status, name, AT_);
}

/**
 * \brief Retrieve an attribute from file (float version).
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
void ParallelIoHdf5::b_getAttribute(T* value, const MString& name, const MString& datasetName,
                                    const size_type totalCount) {
  TRACE();

  // Open the attribute
  hid_t attribute_id;
  if(datasetName == "") { // Attribute is attached to the file
    attribute_id = H5Aopen(b_h5Id, name.c_str(), H5P_DEFAULT);
    b_error(attribute_id, name, AT_);
  } else { // Attribute is attached to a dataset
    attribute_id = H5Aopen_by_name(b_h5Id, datasetName.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    b_error(attribute_id, name, AT_);
  }

  hid_t attribute_type_id = H5Aget_type(attribute_id);
  b_error(attribute_type_id, name, AT_);

  hsize_t length = H5Aget_storage_size(attribute_id) / H5Tget_size(attribute_type_id);
  if(length != (hsize_t)totalCount) {
    TERMM(1, "Requested attribute (" + name
                 + ") has different number of elements than given. Use getAttributeCount() to "
                   "query number of elements first");
  }

  herr_t status = H5Aread(attribute_id, attribute_type_id, value);
  b_error(status, name, AT_);

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Tclose(attribute_type_id);
  b_error(status, name, AT_);
  status = H5Aclose(attribute_id);
  b_error(status, name, AT_);
}


/**
 * \brief Check if object exists. <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[out] attribute exists or not
 * \param[in] name of the attribute
 * \param[in] preceding path of attribute
 *
 * \details
 */
MBool ParallelIoHdf5::b_hasAttribute(const MString& name, const MString& path) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);
  // Handle whether or not path was given
  hid_t loc_id = 0;
  MBool checkGroup = false;
  if(path.empty()) {
    loc_id = b_h5Id;
  } else {
    herr_t pathExists = H5Lexists(b_h5Id, path.c_str(), link_id);
    if(pathExists != 0) {
      // Group path exists, opening it
      loc_id = H5Oopen(b_h5Id, path.c_str(), H5P_DEFAULT);
      b_error(loc_id, path, AT_);
    } else {
      // Group path doesn't exist, then object cannot exist
      H5Pclose(link_id);
      return 0;
    }

    checkGroup = true;
  }

  // Check if object exists
  herr_t status = H5Lexists(loc_id, name.c_str(), link_id);
  if(status == 0) { // Object does not exist, check if attribute exists
    status = H5Aexists(loc_id, name.c_str());
  }

  H5Pclose(link_id);
  if(checkGroup) H5Oclose(loc_id);
  return status;
}

MBool ParallelIoHdf5::b_hasObject(const MString& path) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, path, AT_);
  // Handle whether or not path was given
  MBool exists = false;
  if(path.empty()) {
    exists = true;
  } else {
    herr_t pathExists = H5Lexists(b_h5Id, path.c_str(), link_id);
    if(pathExists != 0) {
      exists = true;
    } else {
      exists = false;
    }
  }

  H5Pclose(link_id);
  return exists;
}


void ParallelIoHdf5::b_getAttributeCount(const MString& name, size_type* totalCount, const MString& datasetName) {
  TRACE();

  // Open the attribute
  hid_t attribute_id;
  if(datasetName.empty()) { // Attribute is attached to the file
    attribute_id = H5Aopen(b_h5Id, name.c_str(), H5P_DEFAULT);
    b_error(attribute_id, name, AT_);
  } else { // Attribute is attached to a dataset
    attribute_id = H5Aopen_by_name(b_h5Id, datasetName.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    b_error(attribute_id, name, AT_);
  }

  hid_t attribute_type_id = H5Aget_type(attribute_id);
  b_error(attribute_type_id, name, AT_);

  hsize_t length = H5Aget_storage_size(attribute_id) / H5Tget_size(attribute_type_id);

  *totalCount = (size_type)length;

  // Close all identifiers which have been created (except b_h5Id)
  herr_t status = H5Tclose(attribute_type_id);
  b_error(status, name, AT_);
  status = H5Aclose(attribute_id);
  b_error(status, name, AT_);
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
MInt ParallelIoHdf5::b_getAttributeType(const MString& name, const MString& datasetName) {
  TRACE();

  // TERMM(1, "untested Hdf5 specific I/O method, please see comment for how"
  //         " to proceed");
  // WARNING: This method is untested, yet. It is implemented due to be
  // consistent with parallel NetCDF. However, it was not used by any of the
  // testcases.If your code uses this part of the code, please make sure that
  // the I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().
  // Determine variable id
  hid_t varId;
  if(datasetName.empty()) { // File attribute
    varId = b_h5Id;
  } else { // Dataset attribute
    varId = H5Dopen2(b_h5Id, datasetName.c_str(), H5P_DEFAULT);
    b_error(varId, name, AT_);
  }

  // Get attribute type
  hid_t attribute_id = H5Aopen(varId, name.c_str(), H5P_DEFAULT);
  b_error(attribute_id, name, AT_);
  hid_t attrType_id = H5Aget_type(attribute_id);
  b_error(attrType_id, name, AT_);

  // Translate native HDF5 data type to ParallelIo type
  MInt typeId;
  if(H5Tequal(attrType_id, H5T_NATIVE_INT) > 0) {
    typeId = PIO_INT;
  } else if(H5Tequal(attrType_id, H5T_NATIVE_DOUBLE) > 0) {
    typeId = PIO_FLOAT;
  } else if(H5Tequal(attrType_id, H5T_NATIVE_LONG) > 0) {
    typeId = PIO_LONG;
  } else if(H5Tequal(attrType_id, H5T_C_S1) > 0) {
    typeId = PIO_STRING;
  } else if(H5Tequal(attrType_id, H5T_NATIVE_UCHAR) > 0) {
    typeId = PIO_UCHAR;
  } else if(H5Tequal(attrType_id, H5T_NATIVE_ULLONG) > 0) {
    typeId = PIO_ULONGLONG;
  } else {
    typeId = PIO_UNKNOWN_TYPE;
  }

  // Close all previously created identifiers (except b_h5Id)
  herr_t status;
  if(varId != b_h5Id) {
    status = H5Dclose(varId);
    b_error(status, name, AT_);
  }
  status = H5Aclose(attribute_id);
  b_error(status, name, AT_);
  status = H5Tclose(attrType_id);
  b_error(status, name, AT_);

  return typeId;
}

//------------------------------------------------------------------------------
// Data mode methods
//------------------------------------------------------------------------------
/**
 * \brief Writes array data to file (String version). <b>[MPI]</b>
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
void ParallelIoHdf5::b_writeArray(const MString* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                  MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                  MPI_Offset diskStride) {
  TRACE();
  TERMM(1, "untested Hdf5 specific I/O method, please see comment for how"
           " to proceed");
  // WARNING: This method is untested, yet. It is implemented due to be
  // consistent with parallel NetCDF. However, it was not used by any of the
  // testcases.If your code uses this part of the code, please make sure that
  // the I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().

  // Get variable id
  hid_t dset_id = H5Dopen2(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  hid_t dspace_id = H5Dget_space(dset_id);
  b_error(dspace_id, name, AT_);
  hid_t dtype_id = H5Tcopy(H5T_C_S1);
  b_error(dtype_id, name, AT_);
  herr_t status = H5Tset_size(dtype_id, H5T_VARIABLE);
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
    // data = tmpScratch.begin();
    data = &tmpScratch[0];
  }

  // Create buffer for writing
  size_type nCount = totalCount * strLen;
  ScratchSpace<char> buf(nCount, FUN_, "buf");

  // Copy data to buffer
  for(MPI_Offset i = 0; i < totalCount; i++) {
    data[i].copy(&buf[i * strLen], strLen, 0);
  }

  hid_t memory_dspace_id = -1;
  // Write array
  if(noChunks == 1) {
    // If number of chunks is one, write everything at once
    memory_dspace_id = H5Screate_simple(noDims, (hsize_t*)count, nullptr);
    b_error(memory_dspace_id, name, AT_);
    if(diskStride == 1) {
      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, nullptr, (hsize_t*)count, nullptr);
    } else {
      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, (hsize_t*)diskStride, (hsize_t*)count,
                                   nullptr);
    }
    b_error(status, name, AT_);
    status = H5Dwrite(dset_id, dtype_id, memory_dspace_id, dspace_id, b_h5DatasetXferHandle, &buf[0]);
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

    ScratchSpace<hsize_t> start_(noDims, FUN_, "start_");
    ScratchSpace<hsize_t> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = max(min(chunkSize, count[0] - i * chunkSize), 0ll);
      const char* buf_ = &buf[0] + min(i * chunkSize * nDSize, (count[0] - 1) * nDSize);
      memory_dspace_id = H5Screate_simple(noDims, &count_[0], nullptr);
      b_error(memory_dspace_id, name, AT_);
      if(diskStride == 1) {
        status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0], nullptr, &count_[0], nullptr);
      } else {
        status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0], (hsize_t*)diskStride, &count_[0], nullptr);
      }
      b_error(status, name, AT_);
      status = H5Dwrite(dset_id, dtype_id, memory_dspace_id, dspace_id, b_h5DatasetXferHandle, &buf_[0]);
      b_error(status, name, AT_);
    }
  }
  // Close all identifiers which have been created (except b_h5Id)
  status = H5Sclose(memory_dspace_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
  status = H5Tclose(dtype_id);
  b_error(status, name, AT_);
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
void ParallelIoHdf5::b_readArray(MString* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                 MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                 MPI_Offset diskStride) {
  TRACE();

  // Get variable id
  hid_t dset_id = H5Dopen(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  hid_t dspace_id = H5Dget_space(dset_id);
  b_error(dspace_id, name, AT_);
  hid_t dtype_id = H5Dget_type(dset_id); // H5Tcopy(H5T_C_S1);


  int length = (int)H5Tget_size(dtype_id);

  if(H5Tget_cset(dtype_id)) {
    dtype_id = H5Tcopy(H5T_C_S1);
    H5Tset_cset(dtype_id, H5T_CSET_UTF8);
  } else {
    dtype_id = H5Tcopy(H5T_C_S1);
  }

  b_error(dtype_id, name, AT_);
  // herr_t status = H5Tset_size(dtype_id, 1);
  herr_t status = H5Tset_size(dtype_id, length);
  b_error(status, name, AT_);

  // Determine total data count
  size_type totalCount = 1;
  for(size_type d = 0; d < noDims - 1; d++) {
    totalCount *= count[d];
  }

  // Determine length of one string
  size_type strLen = length; // count[noDims - 1];

  // Create temporary storage space if needed and set data pointers
  MInt tmpScratchSize = (memoryStride == 1 ? 1 : totalCount);
  ScratchSpace<MString> tmpScratch(tmpScratchSize, FUN_, "tmpStorage");
  MString* data = 0;
  if(memoryStride == 1) {
    data = array;
  } else {
    // data = tmpScratch.begin();
    data = &tmpScratch[0];
  }

  // Create buffer for reading
  size_type nCount = totalCount * strLen;
  ScratchSpace<char> buf(nCount, FUN_, "buf");

  hid_t memory_dspace_id = -1;
  // Read array
  if(noChunks == 1) {
    // If number of chunks is one, write everything at once

    count[0] = 1;

    memory_dspace_id = H5Screate_simple(noDims, (hsize_t*)count, nullptr);
    b_error(memory_dspace_id, name, AT_);
    //    if ( diskStride == 1 ) {
    status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, nullptr, (hsize_t*)count, nullptr);
    //    }
    //    else {
    //      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, (hsize_t*)diskStride,
    //      (hsize_t*)count, nullptr);
    //    }
    b_error(status, name, AT_);
    status = H5Dread(dset_id, dtype_id, memory_dspace_id, dspace_id, H5P_DEFAULT, &buf[0]);
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

    // Determine number of entries for a fixed first dimension index
    size_type nDSize = 1;
    for(size_type d = 1; d < noDims; d++) {
      nDSize *= count[d];
    }

    ScratchSpace<hsize_t> start_(noDims, FUN_, "start_");
    ScratchSpace<hsize_t> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = max(min(chunkSize, count[0] - i * chunkSize), 0ll);
      memory_dspace_id = H5Screate_simple(noDims, &count_[0], nullptr);
      b_error(memory_dspace_id, name, AT_);
      //      if ( diskStride == 1 ) {
      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0], nullptr, &count_[0], nullptr);
      //      }
      //      else {
      //        status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0],
      //                                     (hsize_t*)diskStride, &count_[0], nullptr);
      //      }
      b_error(status, name, AT_);
      status = H5Dread(dset_id, dtype_id, memory_dspace_id, dspace_id, H5P_DEFAULT, &buf[0]);
      b_error(status, name, AT_);

      // Extract strings from buffer
      for(size_type j = start_[0]; j < totalCount; j++) {
        MString tmp;
        tmp.append(&buf[(j - start_[0]) * strLen], strLen);
        data[j].append(tmp.c_str(), 0, strLen);
      }
    }
  }

  // Unpack strided data if necessary
  if(memoryStride != 1) {
    for(MPI_Offset i = 0; i < totalCount; i++) {
      array[memoryStride * i] = tmpScratch[i];
    }
  }

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Sclose(memory_dspace_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
  status = H5Tclose(dtype_id);
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
void ParallelIoHdf5::b_readScalar(MString* scalar, const MString& name) {
  TRACE();

  TERMM(1, "untested Hdf5 specific I/O method, please see comment for how"
           " to proceed");
  // WARNING: This method is untested, yet. It is implemented due to be
  // consistent with parallel NetCDF. However, it was not used by any of the
  // testcases.If your code uses this part of the code, please make sure that
  // the I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().
  // Get variable id
  hid_t dset_id = H5Dopen(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);

  // Read scalar
  hid_t dtype_id = H5Tcopy(H5T_C_S1);
  b_error(dtype_id, name, AT_);
  herr_t status = H5Tset_size(dtype_id, H5T_VARIABLE);
  b_error(status, name, AT_);
  status = H5Dread(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)(scalar));
  b_error(status, name, AT_);

  // Close all identifiers
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
  status = H5Tclose(dtype_id);
  b_error(status, name, AT_);
}


//------------------------------------------------------------------------------
// Attribute methods
//------------------------------------------------------------------------------


/**
 * \brief Creates an attribute in the file (generic version).
 *
 * \author Konstantin Froehlich
 * \date 2015-07-01
 *
 * \param[in] value Attribute value.
 * \param[in] name Attribute name.
 * \param[in] datasetName Dataset (scalar or array) to which the attribute is
 *                        attached (leave empty for file attributes).
 */
template <class T>
void ParallelIoHdf5::b_createAttribute(const T* value, const MString& name, const MString& datasetName, hid_t dtype_id,
                                       const size_type totalCount) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);
  // Allocate dataspace for a scalar
  const auto tmpCnt = (hsize_t)totalCount;
  hid_t dspace_id = (totalCount > 1) ? H5Screate_simple(1, &tmpCnt, nullptr) : H5Screate(H5S_SCALAR);
  b_error(dspace_id, name, AT_);

  hid_t attribute_id = -1;
  herr_t status = -1;

  hid_t loc_id = 0;
  if(!datasetName.empty()) {
    // Create the group if not existent
    status = H5Lexists(b_h5Id, datasetName.c_str(), link_id);
    if(status == 0) {
      hid_t group_id = H5Gcreate(b_h5Id, datasetName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      b_error(group_id, datasetName, AT_);
      H5Gclose(group_id);
    }
    // Open path
    loc_id = H5Oopen(b_h5Id, datasetName.c_str(), H5P_DEFAULT);
    b_error(loc_id, datasetName, AT_);
  } else {
    loc_id = b_h5Id;
  }

  // If the attribute already exists, delete it.
  if(b_hasAttribute(name, datasetName)) {
    if(datasetName.empty()) { // The attribute is attached to the file
      status = H5Adelete(loc_id, name.c_str());
      b_error(status, name, AT_);
    } else { // The attribute is attached to a group
      status = H5Adelete_by_name(loc_id, datasetName.c_str(), name.c_str(), H5P_DEFAULT);
      b_error(status, datasetName, AT_);
    }
  }
  // Create a new attribute
  attribute_id = H5Acreate(loc_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  b_error(attribute_id, name, AT_);

  // Write out the attribute
  status = H5Awrite(attribute_id, dtype_id, value);
  b_error(status, name, AT_);

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Aclose(attribute_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  if(!datasetName.empty()) {
    status = H5Oclose(loc_id);
    b_error(status, name, AT_);
  }
}


/**
 * \brief Creates an attribute in the file (string version).
 *
 * \author Konstantin Froehlich, Felix Wietbuescher
 * \date 2015-07-01
 *
 * \param[in] value Attribute value.
 * \param[in] name Attribute name.
 * \param[in] datasetName Dataset (scalar or array) to which the attribute is
 *                        attached (leave empty for file attributes).
 */
template <>
void ParallelIoHdf5::b_createAttribute(const MString* value, const MString& name, const MString& datasetName,
                                       hid_t dtype_id, const size_type totalCount) {
  TRACE();

  // Allocate dataspace for a scalar
  if(totalCount > 1) {
    mTerm(1, AT_, "Array of strings attributes not yet supported.");
  }

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);

  hsize_t dims = 1;
  hid_t dspace_id = H5Screate_simple(1, &dims, nullptr);
  dtype_id = H5Tcopy(H5T_C_S1);
  MInt length = strlen(value->c_str());
  H5Tset_size(dtype_id, length + 1);
  b_error(dspace_id, name, AT_);

  hid_t attribute_id = -1;
  herr_t status = -1;

  hid_t loc_id = 0;
  if(!datasetName.empty()) {
    // Create the group if not existent
    status = H5Lexists(b_h5Id, datasetName.c_str(), link_id);
    if(status == 0) {
      hid_t group_id = H5Gcreate(b_h5Id, datasetName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      b_error(group_id, datasetName, AT_);
      H5Gclose(group_id);
    }
    // Open path
    loc_id = H5Oopen(b_h5Id, datasetName.c_str(), H5P_DEFAULT);
    b_error(loc_id, datasetName, AT_);
  } else {
    loc_id = b_h5Id;
  }

  // If the attribute already exists, delete it.
  if(b_hasAttribute(name, datasetName)) {
    if(datasetName.empty()) { // The attribute is attached to the file
      status = H5Adelete(loc_id, name.c_str());
      b_error(status, name, AT_);
    } else { // The attribute is attached to a group
      status = H5Adelete_by_name(loc_id, datasetName.c_str(), name.c_str(), H5P_DEFAULT);
      b_error(status, datasetName, AT_);
    }
  }
  // Create a new attribute
  attribute_id = H5Acreate(loc_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  b_error(attribute_id, name, AT_);

  // Write out the attribute
  status = H5Awrite(attribute_id, dtype_id, value->c_str());
  b_error(status, name, AT_);

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Aclose(attribute_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Tclose(dtype_id);
  b_error(status, name, AT_);
  if(!datasetName.empty()) {
    status = H5Oclose(loc_id);
    b_error(status, name, AT_);
  }
}


/**
 * \brief Set an attribute in the file (generic version).
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar
 *
 * \param[in] value Attribute value.
 * \param[in] name Attribute name.
 * \param[in] datasetName Dataset (scalar or array) to which the attribute is
 *                        attached (leave empty for file attributes).
 */
template <class T>
void ParallelIoHdf5::b_setAttribute(const T* value, const MString& name, const MString& datasetName,
                                    const size_type totalCount) {
  TRACE();

  // Create the type and specify its length
  hid_t dtype_id = H5Tcopy(hdf5_traits<T>::type());
  b_error(dtype_id, name, AT_);

  b_createAttribute(value, name, datasetName, dtype_id, totalCount);
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
void ParallelIoHdf5::b_setAttribute(const MString* value, const MString& name, const MString& datasetName,
                                    const size_type totalCount) {
  TRACE();

  if(totalCount > 1) mTerm(1, AT_, "Array of strings attributes not yet supported.");

  // Create the type and specify its length
  hid_t dtype_id = H5Tcopy(H5T_C_S1);
  b_error(dtype_id, name, AT_);
  herr_t status = H5Tset_size(dtype_id, H5T_VARIABLE);
  b_error(status, name, AT_);
  b_createAttribute(value, name, datasetName, dtype_id, 1);
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
void ParallelIoHdf5::b_error(MInt status, const MString& name, const MString& location) {
  TRACE();
  if(status < 0) {
    cerr << endl;
    cerr << "*** ERROR in parallelio_hdf5 ***" << endl;
    cerr << "HDF5 error in function " << location << endl;
    cerr << "HDF5 returns status " << status << endl;
    cerr << "The file/variable/attribute in question was: " << name << endl;
    cerr << endl;
    TERMM(1, "HDF5 error in ParallelIoHdf5.");
  }
}


/**
 * \brief Check the status code of a HDF5 operation and output a meaningful
 *        message (Warning)
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
void ParallelIoHdf5::b_warning(MInt status, const MString& name, const MString& location) {
  TRACE();
  if(status < 0) {
    cerr << endl;
    cerr << "*** HDF5 Warning ***" << endl;
    cerr << "HDF5 warning in function " << location << endl;
    cerr << "HDF5 returns status " << status << endl;
    cerr << "The file/variable/attribute in question was: " << name << endl;
    cerr << endl;
  }
}


//------------------------------------------------------------------------------
// Data mode methods
//------------------------------------------------------------------------------
/**
 * \brief Writes array data to file (Hdf5 version). <b>[MPI]</b>
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
 * \param[in] noChunks Number of Chunks for chunked IO.
 */
template <class T>
void ParallelIoHdf5::b_writeArray(const T* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                  MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                  MPI_Offset diskStride) {
  TRACE();

  // Get variable id
  hid_t dset_id = H5Dopen2(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  hid_t dspace_id = H5Dget_space(dset_id);
  b_error(dspace_id, name, AT_);
  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  // Determine total data count
  size_type totalCount = 1;
  for(size_type d = 0; d < noDims; d++) {
    totalCount *= count[d];
  }

  // Create temporary storage space if needed and set data pointers
  MInt tmpScratchSize = (memoryStride == 1 ? 1 : totalCount);
  ScratchSpace<T> tmpScratch(tmpScratchSize, FUN_, "tmpStorage");

  // Pack strided data
  const T* data = 0;
  if(memoryStride == 1) {
    data = array;
  } else {
    for(MPI_Offset i = 0; i < totalCount; i++) {
      tmpScratch[i] = array[memoryStride * i];
    }
    // data = tmpScratch.begin();
    data = &tmpScratch[0];
  }

  hid_t memory_dspace_id;
  herr_t status = -1;
  // Write array
  if(noChunks == 1) {
    // If number of chunks is one, write everything at once
    memory_dspace_id = H5Screate_simple(noDims, (hsize_t*)count, nullptr);
    b_error(memory_dspace_id, name, AT_);
    if(diskStride == 1) {
      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, nullptr, (hsize_t*)count, nullptr);
    } else {
      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, (hsize_t*)diskStride, (hsize_t*)count,
                                   nullptr);
    }
    b_error(status, name, AT_);
    status = H5Dwrite(dset_id, dtype_id, memory_dspace_id, dspace_id, b_h5DatasetXferHandle, data);
    b_error(status, name, AT_);
  } else {
    TERMM(1, "untested Hdf5 specific I/O method, please see comment for how"
             " to proceed");
    // WARNING: This method is untested, yet. It is implemented due to be
    // consistent with parallel NetCDF. However, it was not used by any of the
    // testcases.If your code uses this part of the code, please make sure that
    // the I/O still works as expected and then remove this warning as well as
    // the subsequent TERMM().
    // Determine variable id
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

    ScratchSpace<hsize_t> start_(noDims, FUN_, "start_");
    ScratchSpace<hsize_t> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = std::min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = std::max(std::min(chunkSize, count[0] - i * chunkSize), 0ll);
      const T* data_ = data + std::min(i * chunkSize * nDSize, (count[0] - 1) * nDSize);
      memory_dspace_id = H5Screate_simple(noDims, &count_[0], nullptr);
      b_error(memory_dspace_id, name, AT_);
      if(diskStride == 1) {
        status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0], nullptr, &count_[0], nullptr);
      } else {
        status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0], (hsize_t*)diskStride, &count_[0], nullptr);
      }
      b_error(status, name, AT_);
      status = H5Dwrite(dset_id, dtype_id, memory_dspace_id, dspace_id, b_h5DatasetXferHandle, data_);
      b_error(status, name, AT_);
    }
  }

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Sclose(memory_dspace_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Tclose(dtype_id);
  b_error(status, name, AT_);
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
}


/**
 * \brief Writes scalar data to file (Hdf5 version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] scalar Value to write.
 * \param[in] name Name of the dataset.
 */
template <class T>
void ParallelIoHdf5::b_writeScalar(T scalar, const MString& name) {
  TRACE();

  // Get variable id
  hid_t dset_id = H5Dopen(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  // Write scalar
  herr_t status = H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scalar);
  b_error(status, name, AT_);

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
}


/**
 * \brief Read array data from file (float version). <b>[MPI]</b>
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
 * \param[in] noChunks Number of Chunks for chunked IO.
 *
 * \details This method is the exact counterpart to h5WriteArray().
 */
template <class T>
void ParallelIoHdf5::b_readArray(T* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                 MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                 MPI_Offset diskStride) {
  TRACE();

  // Get variable id
  hid_t dset_id = H5Dopen(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  hid_t dspace_id = H5Dget_space(dset_id);
  b_error(dspace_id, name, AT_);
  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  // Determine total data count
  size_type totalCount = 1;
  for(size_type d = 0; d < noDims; d++) {
    totalCount *= count[d];
  }

  // Create temporary storage space if needed and set data pointers
  MInt tmpScratchSize = (memoryStride == 1 ? 1 : totalCount);
  ScratchSpace<T> tmpScratch(tmpScratchSize, FUN_, "tmpStorage");
  T* data = 0;
  if(memoryStride == 1) {
    data = array;
  } else {
    // data = tmpScratch.begin();
    data = &tmpScratch[0];
  }

  hid_t memory_dspace_id;
  herr_t status = -1;
  // Read array
  if(noChunks == 1) {
    // If number of chunks is one, write everything at once
    memory_dspace_id = H5Screate_simple(noDims, (hsize_t*)count, nullptr);
    b_error(memory_dspace_id, name, AT_);
    if(diskStride == 1) {
      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, nullptr, (hsize_t*)count, nullptr);
    } else {
      status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, (hsize_t*)start, (hsize_t*)diskStride, (hsize_t*)count,
                                   nullptr);
    }
    b_error(status, name, AT_);
    status = H5Dread(dset_id, dtype_id, memory_dspace_id, dspace_id, H5P_DEFAULT, data);
    b_error(status, name, AT_);
  } else {
    TERMM(1, "untested Hdf5 specific I/O method, please see comment for how"
             " to proceed");
    // WARNING: This method is untested, yet. It is implemented due to be
    // consistent with parallel NetCDF. However, it was not used by any of the
    // testcases.If your code uses this part of the code, please make sure that
    // the I/O still works as expected and then remove this warning as well as
    // the subsequent TERMM().
    // Determine variable id
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

    ScratchSpace<hsize_t> start_(noDims, FUN_, "start_");
    ScratchSpace<hsize_t> count_(noDims, FUN_, "count_");

    std::copy(start, start + noDims, &start_[0]);
    std::copy(count, count + noDims, &count_[0]);

    for(size_type i = 0; i < noChunks; i++) {
      start_[0] = std::min(start[0] + i * chunkSize * diskStride, count[0] * diskStride - 1);
      count_[0] = std::max(std::min(chunkSize, count[0] - i * chunkSize), 0ll);
      T* data_ = data + std::min(i * chunkSize * nDSize, (count[0] - 1) * nDSize);
      memory_dspace_id = H5Screate_simple(noDims, &count_[0], nullptr);
      b_error(memory_dspace_id, name, AT_);
      if(diskStride == 1) {
        status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0], nullptr, &count_[0], nullptr);
      } else {
        status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &start_[0], (hsize_t*)diskStride, &count_[0], nullptr);
      }
      b_error(status, name, AT_);
      status = H5Dread(dset_id, dtype_id, memory_dspace_id, dspace_id, H5P_DEFAULT, data_);
      b_error(status, name, AT_);
    }
  }

  // Unpack strided data if necessary
  if(memoryStride != 1) {
    for(MPI_Offset i = 0; i < totalCount; i++) {
      array[memoryStride * i] = tmpScratch[i];
    }
  }

  // Close all identifiers which have been created (except b_h5Id)
  status = H5Sclose(memory_dspace_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Tclose(dtype_id);
  b_error(status, name, AT_);
  status = H5Dclose(dset_id);
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
template <class T>
void ParallelIoHdf5::b_readScalar(T* scalar, const MString& name) {
  TRACE();

  // Get variable id
  hid_t dset_id = H5Dopen(b_h5Id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  // Read scalar
  herr_t status = H5Dread(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, scalar);
  b_error(status, name, AT_);

  // Close all identifiers
  status = H5Dclose(dset_id);
  b_error(status, name, AT_);
}

/**
 * \brief Read double data from file. <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-09
 *
 * \param[in] array to read
 * \param[in] preceding path of dataset
 * \param[in] name of dataset
 * \param[in] number of dimensions
 * \param[in] offset in different space directions
 * \param[in] number of cells in different space directions
 * * \details
 */
template <class T>
void ParallelIoHdf5::b_readArray(T* array, const MString& path, const MString& name, const size_type noDims,
                                 const size_type* start, const size_type* count) {
  TRACE();

  // Handle whether or not path was given
  hid_t loc_id = 0;
  MInt checkLoc_id = 0;
  if(path.empty()) {
    loc_id = b_h5Id;
  } else {
    loc_id = H5Oopen(b_h5Id, path.c_str(), H5P_DEFAULT);
    b_error(loc_id, path, AT_);
    checkLoc_id = 1;
  }
  // Get variable id
  hid_t dset_id = H5Oopen(loc_id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);
  // Get file space
  hid_t fspace_id = H5Dget_space(dset_id);
  b_error(fspace_id, name, AT_);

  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  // Prepare data
  ScratchSpace<hsize_t> start_(noDims, FUN_, "start_");
  ScratchSpace<hsize_t> count_(noDims, FUN_, "count_");

  std::copy(start, start + noDims, &start_[0]);
  std::copy(count, count + noDims, &count_[0]);

  herr_t status = -1;
  status = H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, &start_[0], NULL, &count_[0], NULL);
  b_error(status, name, AT_);
  hid_t dspace_id = H5Screate_simple(noDims, &count_[0], NULL);
  b_error(dspace_id, name, AT_);
  status = H5Dread(dset_id, dtype_id, dspace_id, fspace_id, H5P_DEFAULT, array);
  b_error(status, name, AT_);

  // Close all indetifiers
  if(checkLoc_id) {
    status = H5Oclose(loc_id);
    b_error(status, name, AT_);
  }
  status = H5Dclose(dset_id);
  b_error(status, path, AT_);
  status = H5Sclose(fspace_id);
  b_error(status, name, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
}


/**
 * \brief Write data array, ghost cell version. <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-09
 *
 * \param[in] array to write
 * \param[in] preceding path of dataset
 * \param[in] name of dataset
 * \param[in] number of dimensions
 * \param[in] offset in different space directions
 * \param[in] number of cells in different space directions
 * \param[in] number of ghost cells
 *
 * \details
 */
template <class T>
void ParallelIoHdf5::b_writeArray(const T* array, const MString& path, const MString& name, const size_type noDims,
                                  const size_type* start, const size_type* count, const size_type* ghost) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);
  // Handle whether or not path was given
  hid_t loc_id = b_h5Id;
  MInt checkLoc_id = 0;
  herr_t status = -1;
  if(path.empty()) {
    loc_id = b_h5Id;
  } else {
    loc_id = H5Oopen(b_h5Id, path.c_str(), H5P_DEFAULT);
    b_error(loc_id, path, AT_);
    checkLoc_id = 1;
  }

  // Check if dataset exists
  status = H5Lexists(loc_id, name.c_str(), link_id);
  b_error(status, name, AT_);
  // Get variable id
  hid_t dset_id = H5Oopen(loc_id, name.c_str(), H5P_DEFAULT);
  b_error(dset_id, name, AT_);

  hid_t dtype_id = H5Dget_type(dset_id);
  b_error(dtype_id, name, AT_);

  // Set memory information
  ScratchSpace<hsize_t> count_(noDims, FUN_, "count_");
  ScratchSpace<hsize_t> start_(noDims, FUN_, "start_");
  ScratchSpace<hsize_t> ghost_(noDims, FUN_, "ghost_");

  for(size_type i = 0; i < noDims; i++) {
    count_[i] = count[i] + (2 * ghost[i]);
  }
  std::copy(start, start + noDims, &start_[0]);
  std::copy(ghost, ghost + noDims, &ghost_[0]);

  hid_t dspace_id = H5Screate_simple(noDims, &count_[0], nullptr);
  b_error(dspace_id, name, AT_);
  // Reassing count values to count_ array before selecting hyperslab
  std::copy(count, count + noDims, &count_[0]);
  status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &ghost_[0], nullptr, &count_[0], nullptr);
  b_error(status, name, AT_);
  // Obtain space
  hid_t fspace_id = H5Dget_space(dset_id);
  b_error(fspace_id, name, AT_);
  // Select hyperslab
  status = H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, &start_[0], nullptr, &count_[0], nullptr);
  b_error(status, name, AT_);
  status = H5Dwrite(dset_id, dtype_id, dspace_id, fspace_id, b_h5DatasetXferHandle, array);

#if !defined(HOST_Hawk)
  // HDF5 installation on Hawk produces strange warnings/errors, although files are ok
  b_error(status, name, AT_);
#endif

  status = H5Pclose(link_id);
  b_error(status, name, AT_);
  if(checkLoc_id) {
    status = H5Oclose(loc_id);
    b_error(status, name, AT_);
  }
  status = H5Dclose(dset_id);
  b_error(status, path, AT_);
  status = H5Sclose(dspace_id);
  b_error(status, name, AT_);
  status = H5Sclose(fspace_id);
  b_error(status, name, AT_);
}


/**
 * \brief Read attribute. <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[in] type of dataset
 * \param[in] preceding path of dataset
 * \param[in] name of dataset
 * \param[in] number of dimensions
 * \param[in] number of cells in different space directions
 *
 * \details
 */
void ParallelIoHdf5::b_defineArray(maiabd_type dsetType, const MString& path, const MString& name,
                                   const size_type noDims, const size_type* count) {
  TRACE();

  hid_t dtype_id;
  switch(dsetType) {
    case PIO_FLOAT: {
      dtype_id = H5T_NATIVE_DOUBLE;
      break;
    }
    case PIO_INT: {
      dtype_id = H5T_NATIVE_INT;
      break;
    }
    case PIO_LONG: {
      dtype_id = H5T_NATIVE_LONG;
      break;
    }
    case PIO_STRING: {
      TERMM(1, "Not yet implemented");
      break;
    }
    case PIO_UCHAR: {
      dtype_id = H5T_NATIVE_UCHAR;
      break;
    }
    case PIO_ULONGLONG: {
      dtype_id = H5T_NATIVE_ULLONG;
      break;
    }

    default: {
      TERMM(1, "Invalid ParallelIo data type!");
    }
  }

  ScratchSpace<hsize_t> datcount(noDims, FUN_, "count_");
  std::copy(count, count + noDims, &datcount[0]);

  hid_t datspace = H5Screate_simple(noDims, &datcount[0], NULL);
  b_error(datspace, name, AT_);
  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);

  // Handle whether or not path was given
  hid_t loc_id = b_h5Id;
  MInt checkLoc_id = 0;
  herr_t status = -1;
  if(path.empty()) {
    loc_id = b_h5Id;
  } else {
    // Create the group if not existent
    status = H5Lexists(b_h5Id, path.c_str(), link_id);
    if(status == 0) {
      // b_createGroup(path);
      loc_id = H5Gcreate(b_h5Id, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      loc_id = H5Gclose(loc_id);
      b_error(loc_id, path, AT_);
    }
    loc_id = H5Oopen(b_h5Id, path.c_str(), H5P_DEFAULT);
    b_error(loc_id, name, AT_);

    checkLoc_id = 1;
  }

  status = H5Lexists(loc_id, name.c_str(), link_id);
  b_error(status, path, AT_);

  hid_t dset_id = H5Dcreate2(loc_id, name.c_str(), dtype_id, datspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  b_error(dset_id, name, AT_);

  status = H5Dclose(dset_id);
  b_error(status, path, AT_);

  // Close all indetifiers
  if(checkLoc_id) {
    status = H5Oclose(loc_id);
    b_error(status, name, AT_);
  }

  status = H5Pclose(link_id);
  b_error(status, name, AT_);
  status = H5Sclose(datspace);
  b_error(status, name, AT_);
}


/**
 * \brief Gets the size of a dataset in space directions. <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[in] name of dataset
 * \param[in] preceding path of dataset
 * \param[in] number of dimensions
 * \param[out] number of cells in different space directions
 *
 * \details
 */
void ParallelIoHdf5::b_getDatasetSize(const MString& name, const MString& path, size_type noDims, size_type* data) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);
  // Handle whether or not path was given
  hid_t loc_id = 0;
  MInt checkLoc_id = 0;
  if(path.empty()) {
    loc_id = b_h5Id;
  } else {
    loc_id = H5Oopen(b_h5Id, path.c_str(), H5P_DEFAULT);
    b_error(loc_id, path, AT_);
    checkLoc_id = 1;
  }

  // Check if object exists
  herr_t status = H5Lexists(loc_id, name.c_str(), link_id);
  if(status == 0) { // Object does not exsit
    m_log << "WARNING: Datatset " << name << " does not exist at path " << path << std::endl;
  } else {
    hid_t data_id = H5Oopen(loc_id, name.c_str(), H5P_DEFAULT);
    // 1) Get space information
    hid_t filespace = H5Dget_space(data_id);
    ScratchSpace<hsize_t> count(noDims, FUN_, "count");
    ScratchSpace<hsize_t> maxcount(noDims, FUN_, "maxcount");
    H5Sget_simple_extent_dims(filespace, &count[0], &maxcount[0]);
    // 2) Get the size
    for(size_type i = 0; i < noDims; i++) {
      data[i] = count[i];
    }
    H5Sclose(filespace);
    H5Oclose(data_id);
  }
  H5Pclose(link_id);
  if(checkLoc_id) H5Oclose(loc_id);
}


/**
 * \brief Returns number of dimensions of dataset <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[out] number of space dimensions of dataset
 * \param[in] name of dataset
 * \param[in] preceding path of dataset
 *
 * \details
 */
ParallelIo::size_type ParallelIoHdf5::b_getDatasetNoDims(const MString& name, const MString& path) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, name, AT_);
  // Handle whether or not path was given
  hid_t loc_id = 0;
  MBool checkLoc_id = false;
  if(path.empty()) {
    loc_id = b_h5Id;
  } else {
    loc_id = H5Oopen(b_h5Id, path.c_str(), H5P_DEFAULT);
    b_error(loc_id, path, AT_);
    checkLoc_id = true;
  }

  hid_t data_id = H5Oopen(loc_id, name.c_str(), H5P_DEFAULT);
  b_error(data_id, name, AT_);

  // 1) Get space information
  hid_t dspace_id = H5Dget_space(data_id);
  b_error(dspace_id, name, AT_);

  size_type noDims = H5Sget_simple_extent_ndims(dspace_id);
  b_error(noDims, name, AT_);

  herr_t status = H5Sclose(dspace_id);
  b_error(status, name, AT_);

  status = H5Oclose(data_id);
  b_error(status, name, AT_);

  status = H5Pclose(link_id);
  b_error(status, name, AT_);

  if(checkLoc_id) {
    status = H5Oclose(loc_id);
    b_error(status, name, AT_);
  }

  return noDims;
}


/**
 * \brief Creates an HDF5 group <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[in] path of group
 *
 * \details
 */
void ParallelIoHdf5::b_createGroup(const MString& path) {
  TRACE();

  hid_t link_id = H5Pcreate(H5P_LINK_ACCESS);
  b_error(link_id, path, AT_);

  // Get path tokens
  std::vector<MString> p;
  tokenize(path, p, "/", false);

  MString parentpath = "/";
  MString objectpath = "/";
  // Get the count number of groups, datasets etc in the path
  MInt num = distance(p.begin(), p.end()) - 1;

  for(MInt i = 0; i < num; ++i) {
    if(p[i] == "") {
      continue;
    }

    if(i != 0) {
      objectpath.append("/");
      parentpath.append(p[i - 1]);
      parentpath.append("/");
    }
    objectpath.append(p[i]);
    herr_t status = H5Lexists(b_h5Id, objectpath.c_str(), link_id);
    b_error(status, path, AT_);

    if(status == 0) { // Path does not exist
      if(i == 0) {
        hid_t group_id = H5Gcreate(b_h5Id, p[i].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        b_error(group_id, path, AT_);
        status = H5Gclose(group_id);
        b_error(status, path, AT_);
      } else { // We are in another group
        hid_t group_id = H5Gopen(b_h5Id, parentpath.c_str(), H5P_DEFAULT);
        b_error(group_id, path, AT_);
        hid_t new_group_id = H5Gcreate(group_id, p[i].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        b_error(new_group_id, path, AT_);

        // Close current groups
        status = H5Gclose(new_group_id);
        b_error(status, path, AT_);
        status = H5Gclose(group_id);
        b_error(status, path, AT_);
      }
    }
  }

  herr_t status = H5Pclose(link_id);
  b_error(status, path, AT_);
}


/**
 * \brief Gets the groups in the given group <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[out] name of the groups inside the given group
 * \param[in] path of group
 * \details
 */
void ParallelIoHdf5::b_getGroupNames(vector<MString>& names, const MString& path) {
  TRACE();
  if(b_getNoGroups(path) <= 0) {
    return;
  }

  MString pathStr = "/";
  if(!path.empty()) {
    pathStr = path;
  }

  hid_t group_id = H5Gopen(b_h5Id, pathStr.c_str(), H5P_DEFAULT);
  b_error(group_id, m_fileName, AT_);

  // Iterate through each object in the group.
  hsize_t noObj = 0;
  herr_t status = H5Gget_num_objs(group_id, &noObj);
  for(MInt idx = 0; idx < (signed)noObj; idx++) {
    int obj_type = H5Gget_objtype_by_idx(group_id, (size_t)idx);
    b_error(obj_type, m_fileName, AT_);
    if(obj_type == H5G_GROUP) {
      char GroupName[NC_MAX_NAME];
      status = H5Gget_objname_by_idx(group_id, (hsize_t)idx, GroupName, NC_MAX_NAME);
      b_error(status, m_fileName, AT_);
      names.emplace_back(GroupName);
    }
  }

  status = H5Gclose(group_id);
  b_error(status, m_fileName, AT_);
}

/**
 * \brief Gets the number of groups in the given group <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[out] number of groups inside given group
 * \param[in] path of group
 * \details
 */
ParallelIo::size_type ParallelIoHdf5::b_getNoGroups(const MString& path) {
  TRACE();

  MString pathStr = "/";
  if(!path.empty()) {
    pathStr = path;
  }

  hid_t group_id = H5Gopen(b_h5Id, pathStr.c_str(), H5P_DEFAULT);
  b_error(group_id, m_fileName, AT_);

  // Iterate through each object in the group.
  hsize_t noObj = 0;
  ParallelIo::size_type noGroups = 0;

  herr_t status = H5Gget_num_objs(group_id, &noObj);
  b_error(status, m_fileName, AT_);
  for(MInt idx = 0; idx < (signed)noObj; idx++) {
    int obj_type = H5Gget_objtype_by_idx(group_id, (size_t)idx);
    b_error(obj_type, m_fileName, AT_);
    if(obj_type == H5G_GROUP) {
      noGroups++;
    }
  }

  status = H5Gclose(group_id);
  b_error(status, m_fileName, AT_);

  return noGroups;
}


/**
 * \brief Gets the dataset names in the given group <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[out] name of datasets in given group
 * \param[in] path of group
 * \details
 */
void ParallelIoHdf5::b_getDatasetNames(vector<MString>& names, const MString& path) {
  TRACE();
  if(b_getNoDatasets(path) <= 0) {
    return;
  }

  MString pathStr = "/";
  if(!path.empty()) {
    pathStr = path;
  }

  hid_t group_id = H5Gopen(b_h5Id, pathStr.c_str(), H5P_DEFAULT);
  b_error(group_id, m_fileName, AT_);

  // Iterate through each object in the group.
  hsize_t noObj = 0;
  herr_t status = H5Gget_num_objs(group_id, &noObj);
  for(MInt idx = 0; idx < (signed)noObj; idx++) {
    int obj_type = H5Gget_objtype_by_idx(group_id, (size_t)idx);
    b_error(obj_type, m_fileName, AT_);
    if(obj_type == H5G_DATASET) {
      char DatasetName[NC_MAX_NAME];
      status = H5Gget_objname_by_idx(group_id, (hsize_t)idx, DatasetName, NC_MAX_NAME);
      b_error(status, m_fileName, AT_);
      names.emplace_back(DatasetName);
    }
  }

  status = H5Gclose(group_id);
  b_error(status, m_fileName, AT_);
}

/**
 * \brief Gets the number of datasetes in the given group <b>[MPI]</b>
 *
 * \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de> (HIWI),
 *         Marian Albers (marian) <marian.albers@aia.rwth-aachen.de>
 * \date 2021-09-18
 *
 * \param[out] number of datasets in given group
 * \param[in] path of group
 * \details
 */
ParallelIo::size_type ParallelIoHdf5::b_getNoDatasets(const MString& path) {
  TRACE();

  MString pathStr = "/";
  if(!path.empty()) {
    pathStr = path;
  }

  hid_t group_id = H5Gopen(b_h5Id, pathStr.c_str(), H5P_DEFAULT);
  b_error(group_id, m_fileName, AT_);

  // Iterate through each object in the group.
  hsize_t noObj = 0;
  ParallelIo::size_type noDatasets = 0;

  herr_t status = H5Gget_num_objs(group_id, &noObj);
  b_error(status, m_fileName, AT_);
  for(MInt idx = 0; idx < (signed)noObj; idx++) {
    int obj_type = H5Gget_objtype_by_idx(group_id, (size_t)idx);
    b_error(obj_type, m_fileName, AT_);
    if(obj_type == H5G_DATASET) {
      noDatasets++;
    }
  }

  status = H5Gclose(group_id);
  b_error(status, m_fileName, AT_);

  return noDatasets;
}

// Explicit instantiations for all supported types
template void ParallelIoHdf5::b_setAttribute(const MFloat* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_setAttribute(const MInt* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_setAttribute(const MLong* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_setAttribute(const MUchar* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_setAttribute(const MUlong* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);

// Explicit instantiations for all supported types
template void ParallelIoHdf5::b_getAttribute(MFloat* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_getAttribute(MInt* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_getAttribute(MLong* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_getAttribute(MUchar* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_getAttribute(MChar* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);
template void ParallelIoHdf5::b_getAttribute(MUlong* value, const MString& name, const MString& datasetName,
                                             const size_type totalCount);

template void ParallelIoHdf5::b_writeArray(const MFloat* array, const MString& name, const size_type noDims,
                                           MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                           const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoHdf5::b_writeArray(const MInt* array, const MString& name, const size_type noDims,
                                           MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                           const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoHdf5::b_writeArray(const MLong* array, const MString& name, const size_type noDims,
                                           MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                           const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoHdf5::b_writeArray(const MUchar* array, const MString& name, const size_type noDims,
                                           MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                           const size_type noChunks, MPI_Offset diskStride);
template void ParallelIoHdf5::b_writeArray(const MUlong* array, const MString& name, const size_type noDims,
                                           MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
                                           const size_type noChunks, MPI_Offset diskStride);

template void ParallelIoHdf5::b_writeArray(const MFloat* array, const MString& path, const MString& name,
                                           const size_type noDims, const size_type* start, const size_type* count,
                                           const size_type* ghost);
template void ParallelIoHdf5::b_writeArray(const MInt* array, const MString& path, const MString& name,
                                           const size_type noDims, const size_type* start, const size_type* count,
                                           const size_type* ghost);


template void ParallelIoHdf5::b_readArray(MFloat* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                          MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                          MPI_Offset diskStride);
template void ParallelIoHdf5::b_readArray(MInt* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                          MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                          MPI_Offset diskStride);
template void ParallelIoHdf5::b_readArray(MLong* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                          MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                          MPI_Offset diskStride);
template void ParallelIoHdf5::b_readArray(MUchar* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                          MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                          MPI_Offset diskStride);
template void ParallelIoHdf5::b_readArray(MUlong* array, const MString& name, const size_type noDims, MPI_Offset* start,
                                          MPI_Offset* count, MPI_Offset memoryStride, const size_type noChunks,
                                          MPI_Offset diskStride);

template void ParallelIoHdf5::b_readArray(MFloat* array, const MString& path, const MString& name, size_type noDims,
                                          const size_type* start, const size_type* count);
template void ParallelIoHdf5::b_readArray(MInt* array, const MString& path, const MString& name, size_type noDims,
                                          const size_type* start, const size_type* count);

template void ParallelIoHdf5::b_readScalar(MFloat* scalar, const MString& name);
template void ParallelIoHdf5::b_readScalar(MInt* scalar, const MString& name);
template void ParallelIoHdf5::b_readScalar(MLong* scalar, const MString& name);
template void ParallelIoHdf5::b_readScalar(MUchar* scalar, const MString& name);
template void ParallelIoHdf5::b_readScalar(MUlong* scalar, const MString& name);

template void ParallelIoHdf5::b_writeScalar(MFloat scalar, const MString& name);
template void ParallelIoHdf5::b_writeScalar(MInt scalar, const MString& name);
template void ParallelIoHdf5::b_writeScalar(MLong scalar, const MString& name);
template void ParallelIoHdf5::b_writeScalar(MUchar scalar, const MString& name);
template void ParallelIoHdf5::b_writeScalar(MUlong scalar, const MString& name);

//#ifdef MAIA_MS_COMPILER
// template void ParallelIoHdf5::b_getAttribute(double* value, const MString& name, const MString& datasetName,
//                                    const size_type totalCount);
// template void ParallelIoHdf5::b_getAttribute(int* value, const MString& name,
//                                             const MString& datasetName,
//                                             const size_type totalCount);
// template void ParallelIoHdf5::b_getAttribute(__int64* value, const MString& name,
//                                             const MString& datasetName,
//                                             const size_type totalCount);
// template void ParallelIoHdf5::b_getAttribute(unsigned __int64* value, const MString& name,
//                                             const MString& datasetName,
//                                             const size_type totalCount);

// template void ParallelIoHdf5::b_readScalar(__int64* scalar, const MString& name);
// template void ParallelIoHdf5::b_readScalar(int* scalar, const MString& name);
// template void ParallelIoHdf5::b_readScalar(long* scalar, const MString& name);
// template void ParallelIoHdf5::b_readScalar(double* scalar, const MString& name);
//
// template void ParallelIoHdf5::b_writeScalar(long scalar, const MString& name);
// template void ParallelIoHdf5::b_writeScalar(__int64 scalar, const MString& name);
// template void ParallelIoHdf5::b_writeScalar(int scalar, const MString& name);
// template void ParallelIoHdf5::b_writeScalar(double scalar, const MString& name);

// template void ParallelIoHdf5::b_writeArray(const int* array, const MString& name, const size_type noDims,
//                                  MPI_Offset* start, MPI_Offset* count, MPI_Offset memoryStride,
//                                  const size_type noChunks, MPI_Offset diskStride);
//
// template void ParallelIoHdf5::b_writeArray(const __int64* array, const MString& name,
//                                           const size_type noDims, MPI_Offset* start,
//                                           MPI_Offset* count, MPI_Offset memoryStride,
//                                           const size_type noChunks, MPI_Offset diskStride);
//
// template void ParallelIoHdf5::b_writeArray(const double* array, const MString& name,
//                                           const size_type noDims, MPI_Offset* start,
//                                           MPI_Offset* count, MPI_Offset memoryStride,
//                                           const size_type noChunks, MPI_Offset diskStride);
//
// template void ParallelIoHdf5::b_writeArray(const unsigned char* array, const MString& name,
//                                           const size_type noDims, MPI_Offset* start,
//                                           MPI_Offset* count, MPI_Offset memoryStride,
//                                           const size_type noChunks, MPI_Offset diskStride);

// template void ParallelIoHdf5::b_readArray(int* array, const MString& name, const size_type noDims,
//                                          MPI_Offset* start, MPI_Offset* count,
//                                          MPI_Offset memoryStride, const size_type noChunks,
//                                          MPI_Offset diskStride);
//
// template void ParallelIoHdf5::b_readArray(__int64* array, const MString& name,
//                                          const size_type noDims,
//                                          MPI_Offset* start, MPI_Offset* count,
//                                          MPI_Offset memoryStride, const size_type noChunks,
//                                          MPI_Offset diskStride);
//
// template void ParallelIoHdf5::b_readArray(double* array, const MString& name,
//                                          const size_type noDims, MPI_Offset* start,
//                                          MPI_Offset* count, MPI_Offset memoryStride,
//                                          const size_type noChunks, MPI_Offset diskStride);
//
// template void ParallelIoHdf5::b_readArray(unsigned char* array, const MString& name,
//                                          const size_type noDims, MPI_Offset* start,
//                                          MPI_Offset* count, MPI_Offset memoryStride,
//                                          const size_type noChunks, MPI_Offset diskStride);
//#endif

#endif
