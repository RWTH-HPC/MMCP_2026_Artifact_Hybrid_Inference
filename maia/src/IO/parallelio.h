// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PARALLELIO_H
#define PARALLELIO_H

#include <set>
#include <sys/stat.h>
#include <utility>
#include <vector>
#include "COMM/globalmpiinfo.h"
#include "COMM/mpioverride.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "UTIL/timer.h"
#include "config.h"
#include "typetraits.h"

namespace maia {
namespace parallel_io {

// Typedef for constants
// Used for all constants that are used to interact with ParallelIoPNetcdf
using maiabd_type = MInt;

// File mode constants
//!< File mode to create a new file. Aborts if file already exists
const MInt PIO_CREATE = 0;
// File mode to create a new file. If file already exists, it is overwritten
const MInt PIO_REPLACE = 1;
// File mode to open an existing file for reading and writing
const MInt PIO_APPEND = 2;
// File mode to open an existing file for reading only
const MInt PIO_READ = 3;

// Data type constants
// Specifies all unknown (= unsupported) data types.
const MInt PIO_UNKNOWN_TYPE = -1;
// Specifies data of type MFloat.
const MInt PIO_FLOAT = 0;
// Specifies data of type MInt.
const MInt PIO_INT = 1;
// Specifies data of type MLong.
const MInt PIO_LONG = 2;
// Specifies data of type MString.
const MInt PIO_STRING = 3;
// Specifies data of type MChar.
const MInt PIO_UCHAR = 4;
// Specifies data of type MUlong.
const MInt PIO_ULONGLONG = 5;

} // namespace parallel_io
} // namespace maia

/**
 * \brief This class is intended to do all the heavy lifting when it comes to
 *        reading and writing "big data" files. It supports multiple backends,
 *        i.e. different I/O libraries. The default is "Parallel netCDF".
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-03-27
 *
 * As of fall 2015 there are two backends available: Parallel netCDF and HDF5.
 *
 * The main intention of the class is : complete separation of I/O
 * calls from the user, so that nobody has to deal with the different ways of
 * how to use (Parallel) NetCDF, HDF5, MPI I/O etc. The class also offers an
 * easy way to make sure that data files (i.e. restart files, solution files)
 * and grid files always have the correct format.
 *
 * *Note:* All methods that may be slow or incur a significant overhead since
 * they use MPI communication (either for I/O or other purposes) are marked
 * with <b>[MPI]</b>. While they might not *always* have this overhead, it is
 * better to assume they have and to use them only when necessary.
 *
 * *Note:* Almost all methods in this class are collective, i.e. they have to
 * be called on all ranks within the given MPI communicator to work properly
 * (even if they themselves are not marked with <b>[MPI]</b>). Some methods may
 * fail if you try to call them from a subset of the ranks, but do not rely on
 * a meaningful error.
 *
 * *Note:* Since I/O is often an error-prone area, it is usually advisable to
 * check all I/O operations for possible occurring errors. In order to relieve
 * the user from the burden of doing correct and repeated checking of possible
 * error values/message, ParallelIo adheres to the following rule:
 *
 * <b>ALL ERRORS ARE FATAL!</b>
 *
 * Exceptions to this rule are marked as such. This has the advantage that
 * there is no need to do extra checking in the production code of MAIA itself
 * on the one hand, and on the other hand it is assured that it is not possible
 * that I/O operations leave MAIA in an undefined state if an error occurred.
 */
template <class Backend>
class ParallelIoBase {
  // CRTP
 private:
  Backend& derived() { return static_cast<Backend&>(*this); }
  const Backend& derived() const { return static_cast<const Backend&>(*this); }

  // Public types
 public:
  /**
   * \brief Type used for all size- and offset-related values.
   *
   * This type will be used internally and externally for all values that can
   * be used as size values or offsets. Examples are 'count', 'size', 'stride'
   * etc. The usage of this type is to ensure that this type may be changed
   * with little to no necessary changes within the code.  A change might be
   * necessary e.g. if the previous type's range becomes to small to hold all
   * relevant size information.
   *
   * *Note:* If you change this type, you also MUST check if the static member
   * s_mpiSizeType needs to be changed as well, since it is used as the MPI
   * type to communicate values of the 'size_type' kind.
   */
  using size_type = MLong;

  /**
   * \brief Type used for all constants that are defined in maia::parallel_io.
   */
  using maiabd_type = maia::parallel_io::maiabd_type;

  // Public member functions
 public:
  // Static file system-related methods
  static MBool fileExists(const MString& name, const MPI_Comm& mpiComm);
  static void deleteFile(const MString& name);
  static MBool isParallelIoFile(const MString& name, const MPI_Comm& mpiComm);
  static MString fileExt();

  // Constructors & destructor; protected, so that the user can't create it
  // directly. Only objects of the derived classes should exist
 protected:
  ParallelIoBase(MString fileName, MInt fileMode, const MPI_Comm& mpiComm);
  virtual ~ParallelIoBase();

  // The following methods represent the interface used by the derived classes
 public:
  // Methods to get information about file contents
  MBool hasDataset(const MString& name, MInt dimension);
  MBool hasDataset(const MString& name, const MString& path);
  MBool hasDataset(const MString& name);
  MBool hasObject(const MString& path);
  MInt getDatasetType(const MString& name);
  std::vector<MString> getDatasetNames(const size_type dimensions = -1);
  std::vector<MString> getDatasetNames(const MString& path);
  std::vector<MString> getGroupNames(const MString& path);
  size_type getDatasetNoDims(const MString& name);
  size_type getDatasetNoDims(const MString& name, const MString& path);
  std::vector<size_type> getArrayDims(const MString& name);
  size_type getArraySize(const MString& name, const size_type dimensionId = 0);
  void getArraySize(const MString& name, const MString& path, size_type* datasetSize);
  MBool hasAttribute(const MString& name, const MString& path = "");
  maiabd_type getAttributeType(const MString& name, const MString& datasetName = "");

  // Methods to define arrays and scalars
  void defineArray(maiabd_type type, const MString& name, size_type totalCount);
  void defineArray(maiabd_type type, const MString& name, size_type noDims, size_type* totalCount);
  void defineArray(maiabd_type type, const std::vector<MString>& name, size_type totalCount);
  void defineArray(maiabd_type type, const MString& datasetName, const MString& name, size_type noDims,
                   size_type* totalCount);
  void defineScalar(maiabd_type type, const MString& name);

  // Methods to write data or attributes
  template <class T>
  void writeArray(const T* array, const MString& name, size_type memoryStride = -1, size_type diskStride = -1);
  template <class T>
  void writeArray(const T* array, const std::vector<MString>& names, size_type memoryStride = -1,
                  size_type diskStride = -1);
  template <class T>
  void writeArray(const T* array, const MString& path, const MString& name, const size_type noDims, size_type* offset,
                  size_type* localCount);
  template <class T>
  void writeArray(const T* array, const MString& path, const MString& name, const size_type noDims, size_type* offset,
                  size_type* localCount, size_type* ghost);
  template <class T>
  void writeScalar(T scalar, const MString& name);
  template <class T>
  void setAttribute(const T& value, const MString& name, const MString& datasetName = "");
  void setAttribute(const MChar* value, const MString& name, const MString& datasetName = "");
  template <class T>
  void setAttributes(const T* value, const MString& name, size_type totalCount, const MString& datasetName = "");

  // Methods to read data or attributes
  template <class T>
  void readArray(T* array, const MString& name, size_type memoryStride = -1, size_type diskStride = -1);
  template <class T>
  void readArray(std::vector<T>& array, const MString& name, size_type memoryStride = -1, size_type diskStride = -1) {
    readArray(array.data(), name, memoryStride, diskStride);
  }
  template <class T>
  void readArray(T* array, const std::vector<MString>& names, size_type memoryStride = -1, size_type diskStride = -1);

  template <class T>
  void readArray(T* array, const MString& datasetName, const MString& name, const size_type noDims,
                 const size_type* offset, const size_type* localCount);
  template <class T>
  void readArray(T* array, const MString& datasetName, const MString& name, const size_type noDims,
                 const size_type* localCount);

  template <class T>
  void readScalar(T* scalar, const MString& name);
  template <class T>
  void getAttribute(T* value, const MString& name, const MString& datasetName = "");
  template <class T>
  void getAttribute(T* value, const MString& name, const size_type totalCount, const MString& datasetName = "");
  void getAttributeCount(const MString& name, size_type* totalCount, const MString& datasetName = "");

  // Methods to calculate and set offset
  void setOffset(const size_type localCount, const size_type offset, const size_type noDims, const size_type noChunks);
  void setOffset(const size_type localCount, const size_type offset, const size_type noDims);
  void setOffset(size_type localCount, size_type offset);

  static void calcOffset(size_type localCount, size_type* offset, size_type* totalCount, const MPI_Comm& mpiComm,
                         size_type* noChunks = nullptr);

  // Member variables and methods which are not part of the interface
 protected:
  // Auxiliary methods
  void validateType(const maiabd_type type) const;

  // File related
  const MString m_fileName{};
  const MInt m_fileMode = -1;
  MBool m_isHeaderSaved = false;

  // Other core member variables used for information
  // which will be added into the file header
  MInt m_domainId = -1;
  MInt m_noDomains = -1;

  // MPI-related
  MPI_Comm m_mpiComm{};
  MPI_Info m_mpiInfo = MPI_INFO_NULL;

  // Offset array-related (usually: offset size = number of cells)
  MBool m_offsetsAreSet = false;
  std::vector<size_type> m_localCount = {-1};
  std::vector<size_type> m_offset = {-1};
  size_type m_noChunks = -1;

  // Creation time stamp to output total file lifetime
  MFloat m_creationTime = MFloatNaN;

#ifdef MAIA_EXTRA_DEBUG
  // Debugging related
  std::set<MString> m_unwrittenArrays;
  std::set<MString> m_unwrittenScalars;
  std::set<MString> m_writtenArrays;
  std::set<MString> m_writtenScalars;
  std::set<MString> m_readArrays;
  std::set<MString> m_readScalars;
#endif

  // Maximum allowable size for writing to an array at once (limited by MPI I/O,
  // as the 'count' parameter in the relevant MPI call is a 32-bit integer). The
  // value used here is 2^31 - 3, the '-3' being an extra buffer
  const static size_type s_maxChunkSize = 2147483645;

  // Internal pseudo-datatypes used for MPI-related tasks
  // Used to communicate values of type 'size_type'. Therefore the chosen type
  // must always be big enough to hold a single value of 'size_type'.
  const static MPI_Datatype s_mpiSizeType;
};

// For future Developers of ParallelIo:
// Why here? These headers contain definitions of derived classes of
// ParallelIoBase(-> They can't be on the top). After these includes the
// definitions of member functions of
// ParallelIoBase call functions of these derived classes(-> They have to be
// included before member function definitions and after defining the Base
// Class).
#if !defined(MAIA_WINDOWS)
#include "parallelio_pnetcdf.h"
#elif !defined(WITH_HDF5)
#error Cannot compile on Windows if "WITH_HDF5" is not set
#endif
#if defined(WITH_HDF5)
#include "parallelio_hdf5.h"
#endif

// Add convenience type definition to allow simple usage of the default backend
// "configure.py ? ?" -> PARALLELIO_DEFAULT_BACKEND = ParallelIoPNetcdf
// "configure.py ? ? --io hdf5" -> PARALLELIO_DEFAULT_BACKEND = ParallelIoHdf5
using ParallelIo = PARALLELIO_DEFAULT_BACKEND;

// Set static members
template <typename Backend>
const MPI_Datatype
    ParallelIoBase<Backend>::s_mpiSizeType = maia::type_traits<ParallelIoBase<Backend>::size_type>::mpiType();

//------------------------------------------------------------------------------
// Static file system-related methods
//------------------------------------------------------------------------------


/**
 * \brief Check whether a file exists. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] mpiComm The MPI communicator to be used to check the file
 *                    existence.
 *
 * \return True, if the specified file exists.
 *
 * \details This must be a collective call from all ranks in mpiComm.
 */
template <typename Backend>
MBool ParallelIoBase<Backend>::fileExists(const MString& name, const MPI_Comm& mpiComm) {
  TRACE();

  // Try to stat file on rank 0 - will return something other than zero of file
  // does not exist - and broadcast to others
  int rank, status;
  MPI_Comm_rank(mpiComm, &rank);
  if(rank == 0) {
    struct stat buffer {};
    status = static_cast<int>(stat(name.c_str(), &buffer) == 0);
  }
  MPI_Bcast(&status, 1, MPI_INT, 0, mpiComm, AT_, "status");

  return status != 0;
}


/**
 * \brief Deletes the specified file. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-01
 *
 *
 * \details This method does not check if the file exists before attempting to
 *          delete it.
 */
template <typename Backend>
void ParallelIoBase<Backend>::deleteFile(const MString& name) {
  TRACE();

  MPI_File_delete(const_cast<MChar*>(name.c_str()), MPI_INFO_NULL);
}


/**
 * \brief Check if specified file is a valid ParallelIoHdf5 file. <b>[MPI]</b>
 *
 * \author Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-03-01
 *
 * \param[in] mpiComm MPI communicator to be used.
 *
 * \details This must be a collective call from all ranks in mpiComm.
 */
template <class Backend>
MBool ParallelIoBase<Backend>::isParallelIoFile(const MString& name, const MPI_Comm& mpiComm) {
  TRACE();

  return Backend::b_isValidFile(name, mpiComm);
}


/**
 * \brief Returns backend-specific ending of filename (either ".Netcdf" or
 *  ".Hdf5")
 * \author Konstantin Froehlich
 * \date 2015-10-21
 *
 */
template <class Backend>
MString ParallelIoBase<Backend>::fileExt() {
  TRACE();

  return Backend::b_fileExt();
}


//------------------------------------------------------------------------------
// Constructors & destructor
//------------------------------------------------------------------------------


/**
 * \brief Creates a new object to read and write *big* data files. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-03-26
 * \author (Modified by: )Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-03-01
 *
 * \param[in] fileMode The file mode can be either maia::parallel_io::PIO_CREATE,
 *                     maia::parallel_io::PIO_APPEND, or maia::parallel_io::PIO_READ.
 * \param[in] mpiComm The MPI communicator that should be used to open/create
 *                    the file.
 */
template <class Backend>
ParallelIoBase<Backend>::ParallelIoBase(MString fileName, MInt fileMode, const MPI_Comm& mpiComm)
  : m_fileName(std::move(fileName)), m_fileMode(fileMode), m_mpiComm(mpiComm) {
  TRACE();

  using namespace maia::parallel_io;
  // Make sure that the file mode is supported
  if(fileMode != PIO_CREATE && fileMode != PIO_APPEND && fileMode != PIO_REPLACE && fileMode != PIO_READ) {
    TERMM(1, "File mode must be one of PIO_CREATE/PIO_APPEND/PIO_REPLACE/PIO_READ!");
  }

  // Set domain id and number of domains
  MPI_Comm_rank(m_mpiComm, &m_domainId);
  MPI_Comm_size(m_mpiComm, &m_noDomains);

  // Store creation time stamp
  m_creationTime = wallTime();
}


/**
 * \brief empty destructor. <b>[MPI]</b>
 *
 * \author Konstantin Froehlich
 * \date 2013-03-26
 */
template <class Backend>
ParallelIoBase<Backend>::~ParallelIoBase() {
  TRACE();

  // Log the total file lifetime
  const MFloat lifetime = wallTime() - m_creationTime;

  const MInt maxLineLength = 256;
  MChar b[maxLineLength];
  snprintf(b, maxLineLength, "=== PARALLELIO FILE LIFETIME: %-35s | %.4e s |\n", m_fileName.c_str(), lifetime);
  if(m_domainId == 0) {
    std::cerr << b;
  }
  m_log << b;
}


//------------------------------------------------------------------------------
// Methods to get information about file contents
//------------------------------------------------------------------------------


/**
 * \brief Check if the file contains an dataset with the given name and
 * dimension.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name The name of the dataset that should be checked.
 * \param[in] dimension Dimension of dataset to match, -1 will
 * match any dataset, 0 will only match scalars, 1 arrays with one dimensions,
 * ...
 *
 * \return True if an dataset of the given name with given dimension exists.
 */
template <class Backend>
MBool ParallelIoBase<Backend>::hasDataset(const MString& name, const MInt dimension) {
  TRACE();

  MBool returnValue = derived().b_hasDataset(name, dimension);

  return returnValue;
}


/**
 * \brief Check if the file contains an dataset with the given name
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name The name of the dataset that should be checked.
 *
 * \return True if an dataset of the given name exists.
 */
template <class Backend>
MBool ParallelIoBase<Backend>::hasDataset(const MString& name) {
  TRACE();

  MBool returnValue = derived().b_hasDataset(name, -1);

  return returnValue;
}

template <class Backend>
MBool ParallelIoBase<Backend>::hasDataset(const MString& name, const MString& path) {
  TRACE();

  MBool returnValue = derived().b_hasDataset(name, path);

  return returnValue;
}

template <class Backend>
MBool ParallelIoBase<Backend>::hasObject(const MString& name) {
  TRACE();

  MBool returnValue = derived().b_hasObject(name);

  return returnValue;
}


/**
 * \brief Returns the data type of an array.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-01
 * \author Modified by: Ramandeep Jain (HiWi) <ramandeepjain@gmail.com>, Konstantin Froehlich
 * \date 2015-Mar
 *
 * \param[in] name Name of the array to check.
 *
 * \return The data type as defined in the namespace maia::parallel_io.
 *
 * \details Calls hasDataset() first to check if the array exists and aborts
 *          fatally if not.
 */
template <class Backend>
MInt ParallelIoBase<Backend>::getDatasetType(const MString& name) {
  TRACE();

  // Abort if dataset does not exist
  if(!hasDataset(name)) {
    TERMM(1, "The specified dataset '" + name + "' does not exist.");
  }

  return derived().b_getDatasetType(name);
}


/**
 * \brief Returns a vector with the names of all existing datasets with given
 *        dimensionality in the file.
 *
 * \author Ansgar Niemoeller, , Konstantin Froehlich
 * \date 2015-08-19
 *
 * \param[in] dimension Return datasets with given dimensionality (default -1)
 * \param[out] names Contains all found dataset names with the given
 *             dimensionality (if any).
 */
template <class Backend>
std::vector<MString> ParallelIoBase<Backend>::getDatasetNames(const size_type dimension) {
  TRACE();

  std::vector<MString> names;
  derived().b_getDatasetNames(names, dimension);

  return names;
}


template <class Backend>
std::vector<MString> ParallelIoBase<Backend>::getDatasetNames(const MString& path) {
  TRACE();

  std::vector<MString> names;
  derived().b_getDatasetNames(names, path);

  return names;
}

template <class Backend>
std::vector<MString> ParallelIoBase<Backend>::getGroupNames(const MString& path) {
  TRACE();

  std::vector<MString> names;
  derived().b_getGroupNames(names, path);

  return names;
}

/**
 * \brief Get the number of dimensions of a dataset with given name
 *
 * \author Ansgar Niemoeller, Konstantin Froehlich
 * \date 2015-08-19
 *
 * \param[in] name Name of dataset to check
 * \param[out] noDims Number of dimensions of the dataset
 */
template <class Backend>
typename ParallelIoBase<Backend>::size_type ParallelIoBase<Backend>::getDatasetNoDims(const MString& name) {
  TRACE();

  // Abort if dataset does not exist
  if(!hasDataset(name)) {
    TERMM(1, "The specified dataset '" + name + "' does not exist.");
  }

  size_type noDims = derived().b_getDatasetNoDims(name);

  return noDims;
}

template <class Backend>
typename ParallelIoBase<Backend>::size_type ParallelIoBase<Backend>::getDatasetNoDims(const MString& name,
                                                                                      const MString& path) {
  TRACE();

  // Abort if dataset does not exist
  if(!hasDataset(name, path)) {
    TERMM(1, "The specified dataset '" + name + "' does not exist.");
  }

  size_type noDims = derived().b_getDatasetNoDims(name, path);

  return noDims;
}


/**
 * \brief Get the lengths of all dimensions of a dataset with given name
 *
 * \author Ansgar Niemoeller, Konstantin Froehlich
 * \date 2015-08-19
 *
 * \param[in] name Name of dataset to check
 * \param[out] dims Contains the lengths of all dimensions of the dataset
 */
template <class Backend>
std::vector<typename ParallelIoBase<Backend>::size_type> ParallelIoBase<Backend>::getArrayDims(const MString& name) {
  TRACE();

  // Abort if array does not exist
  if(!hasDataset(name) || hasDataset(name, 0)) {
    TERMM(1, "The specified array '" + name + "' does not exist.");
  }

  // Get number of dimensions
  size_type noDims = getDatasetNoDims(name);

  std::vector<size_type> dims;

  // Get lengths of all dimensions
  for(size_type dimId = 0; dimId < noDims; dimId++) {
    dims.push_back(getArraySize(name, dimId));
  }

  return dims;
}


/**
 * \brief Get the length of an array in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-01
 *
 * \param[in] name Name of the array to check.
 * \param[in] dimensionId Dimension of the array to check (default 0).
 *
 * \return Total number of elements of the array in the given dimension.
 */
template <class Backend>
typename ParallelIoBase<Backend>::size_type ParallelIoBase<Backend>::getArraySize(const MString& name,
                                                                                  const size_type dimensionId) {
  TRACE();

  // Abort if dataset does not exist or if its a scalar
  if(!hasDataset(name) || hasDataset(name, 0)) {
    TERMM(1, "The specified array '" + name + "' does not exist.");
  }

  // Abort if dimensionId exceeds number of dimensions
  size_type noDims = getDatasetNoDims(name);
  if(dimensionId >= noDims) {
    TERMM(1, "The specified dimension for array '" + name + "'does not exist.");
  }

  size_type returnValue = derived().b_getDatasetSize(name, dimensionId);

  return returnValue;
}

template <class Backend>
void ParallelIoBase<Backend>::getArraySize(const MString& name, const MString& path, size_type* datasetSize) {
  TRACE();

  // Abort if dataset does not exist or if its a scalar
  if(!hasDataset(name, path)) {
    TERMM(1, "The specified array '" + name + "' does not exist.");
  }

  size_type noDims = derived().b_getDatasetNoDims(name, path);
  derived().b_getDatasetSize(name, path, noDims, datasetSize);
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
template <class Backend>
MBool ParallelIoBase<Backend>::hasAttribute(const MString& name, const MString& path) {
  TRACE();

  // If the dataset name is empty, this means a file attribute should be
  // checked.
  // Otherwise a dataset attribute (scalar or array) is to be checked - abort if
  // no dataset can be found.
  if(!path.empty()) {
    if(!hasObject(path)) {
      return false;
    }
  }

  MBool returnValue = derived().b_hasAttribute(name, path);

  return returnValue;
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
 *
 * \details Calls hasAttribute() first to check if the attribute exists and
 *          aborts fatally if not.
 */
template <class Backend>
MInt ParallelIoBase<Backend>::getAttributeType(const MString& name, const MString& datasetName) {
  TRACE();

  // Abort if attribute does not exist
  if(!hasAttribute(name, datasetName)) {
    mTerm(1, AT_, "The specified attribute does not exist.");
  }

  return derived().b_getAttributeType(name, datasetName);
}


//------------------------------------------------------------------------------
// Methods to define arrays and scalars
//------------------------------------------------------------------------------


/**
 * \brief Create a new array in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-08
 *
 * \param[in] type Data type of the array (either maia::parallel_io::PIO_FLOAT,
 *                 maia::parallel_io::PIO_INT, maia::parallel_io::PIO_LONG or
 *                 maia::parallel_io::PIO_STRING or maia::parallel_io::PIO_UCHAR).
 * \param[in] name Name of the array (may not be empty).
 * \param[in] totalCount Total size of the array (i.e. combined size of all
 *                       involved domains).
 *
 * \details This method may not be called after the first call to writeArray()
 *          or writeScalar().
 */
template <class Backend>
void ParallelIoBase<Backend>::defineArray(maiabd_type type, const MString& name, size_type totalCount) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  defineArray(type, name, 1, &totalCount);
}


/**
 * \brief Create a new array in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-08
 *
 * \param[in] type Data type of the array (either maia::parallel_io::PIO_FLOAT,
 *                 maia::parallel_io::PIO_INT, maia::parallel_io::PIO_LONG or
 *                 maia::parallel_io::PIO_STRING or maia::parallel_io::PIO_UCHAR).
 * \param[in] name Name of the array (may not be empty).
 * \param[in] noDims Number of array dimensions
 * \param[in] totalCount Total size of the array for all dimensions (i.e.
 *                       combined size of all involved domains).
 *
 * \details This method may not be called after the first call to writeArray()
 *          or writeScalar().
 */
template <class Backend>
void ParallelIoBase<Backend>::defineArray(maiabd_type type, const MString& name, size_type noDims,
                                          size_type* totalCount) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  using namespace maia::parallel_io;

  // Prevent adding data in read-only mode
  if(m_fileMode == PIO_READ) {
    mTerm(1, AT_, "Cannot add data in read-only mode!");
  }

  // Abort if header was already written
  if(m_isHeaderSaved) {
    mTerm(1, AT_, "Cannot add new data after header was written!");
  }

#ifdef MAIA_EXTRA_DEBUG
  m_unwrittenArrays.insert(name);
#endif

  // Ensure reasonable parameter values
  validateType(type);
  if(name.empty()) {
    TERMM(1, "Invalid name! Name must not be empty.");
  }
  if(noDims <= 0) {
    TERMM(1, "Invalid number of dimensions! Must at least be 1.");
  }

  derived().b_defineArray(type, name, noDims, totalCount);
}


/**
 * \brief Create multiple new arrays in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-08
 *
 * \param[in] type Data type of the arrays (either maia::parallel_io::PIO_FLOAT or
 *                 maia::parallel_io::PIO_INT).
 * \param[in] names Names of the arrays (may not be empty).
 * \param[in] totalCount Total size of each array (i.e. combined size of all
 *                       involved domains).
 *
 * \details This method may not be called after the first call to writeArray()
 *          or writeScalar().
 */
template <class Backend>
void ParallelIoBase<Backend>::defineArray(maiabd_type type, const std::vector<MString>& names, size_type totalCount) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  for(const auto& name : names) {
    defineArray(type, name, totalCount);
  }
}

/**
 * \brief Creates new array in the file.
 *
 * \author Marian Albers <m.albers@aia.rwth-aachen.de
 * \date 2022-09-25
 *
 * \param[in] type Data type of the arrays (either maia::parallel_io::PIO_FLOAT or
 *                 maia::parallel_io::PIO_INT).
 * \param[in] Name of the array
 * \param[in] Number of dimensions of array
 * \param[in] Size of array in each space direction
 *
 */
template <class Backend>
void ParallelIoBase<Backend>::defineArray(maiabd_type type, const MString& datasetName, const MString& name,
                                          size_type noDims, size_type* totalCount) {
  TRACE();

  using namespace maia::parallel_io;

  // Prevent adding data in read-only mode
  if(m_fileMode == PIO_READ) {
    mTerm(1, AT_, "Cannot add data in read-only mode!");
  }

#ifdef MAIA_EXTRA_DEBUG
  m_unwrittenArrays.insert(name);
#endif

  // Ensure reasonable parameter values
  validateType(type);
  if(name.empty()) {
    TERMM(1, "Invalid name! Name must not be empty.");
  }
  if(noDims <= 0) {
    TERMM(1, "Invalid number of dimensions! Must at least be 1.");
  }

  derived().b_defineArray(type, datasetName, name, noDims, totalCount);
}


/**
 * \brief Create a new scalar in the file.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-08
 *
 * \param[in] type Data type of the scalar (either maia::parallel_io::PIO_FLOAT or
 *                 maia::parallel_io::PIO_INT).
 * \param[in] name Name of the scalar (may not be empty).
 *
 * \details This method may not be called after the first call to writeArray()
 *          or writeScalar().
 */
template <class Backend>
void ParallelIoBase<Backend>::defineScalar(maiabd_type type, const MString& name) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  using namespace maia::parallel_io;

  // Prevent adding data in read-only mode
  if(m_fileMode == PIO_READ) {
    TERMM(1, "Cannot add data in read-only mode!");
  }

  // Abort if header was already written
  if(m_isHeaderSaved) {
    TERMM(1, "Cannot add new data after header was written!");
  }

#ifdef MAIA_EXTRA_DEBUG
  m_unwrittenScalars.insert(name);
#endif

  // Ensure reasonable parameter values
  validateType(type);
  if(name.empty()) {
    TERMM(1, "Invalid name! Name must not be empty.");
  }

  derived().b_defineScalar(type, name);
}


//------------------------------------------------------------------------------
// Methods to write data or attributes
//------------------------------------------------------------------------------


/**
 * \brief Write array data to file. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the array (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data in memory.
 * \param[in] name Name of the array in the file.
 * \param[in] memoryStride Offset between two data values in *memory* (or -1 to set
 *                   automatically).
 * \param[in] diskStride Offset between two data values on *disk* (or -1 to set
 *                   automatically).
 *
 *  The size and file offset have to be set using setOffset() before this
 *  method is called.  If no stride is given (or -1), the stride will be set to
 *  1.
 *
 *  Example:
 *
 *  Suppose you would like to store all z-coordinates of your cells in a file,
 *  and assign them the name 'coordinates_2'. You have a total of 1000 cells
 *  which are split across two domains, domain0 with 550 cells and domain1 with
 *  the remaing 450.
 *
 *  Then you first have to define the array with a call
 *
 *      defineArray(PIO_FLOAT, "coordinates_2", 1000);
 *
 *  Note that '1000' denoting the *total* count, i.e. the number of values you
 *  would like to write from all domains (or, alternatively worded, this is the
 *  size the dataset will have in the file). If you do not know the total
 *  count yet, you can use calcOffset() to conveniently calculate it.
 *
 *  Next you have to set the offset, so each processor knows where it should
 *  start writing.  In this case its rather simple: domain0 should write its
 *  first value at the beginning of the file, and domain1 should leave a gap of
 *  550 values for domain0. Thus you have to call setOffset() once on each
 *  domain, i.e.
 *
 *      setOffset(550, 0);   // <-- This line on domain0. We have 550 cells and
 *                           //     start at the beginning.
 *      setOffset(450, 550); // <-- This line on domain1. We have 450 cells and
 *                           //     leave a gap of 550.
 *
 *  Now you can start writing. For this you need to know how your data is
 *  distributed in memory. In MAIA, most cell data is stored contiguously in
 *  memory, grouped by cell. This means that for 3D, the coordinates are layed
 *  out in memory as (noted as cellId:direction)
 *
 *      | 0:x | 0:y | 0:z | 1:x | 1:y | 1:z | 2:x | 2:y | 2:z | 3:x | 3:y | ...
 *
 *  For the writeArray() methods you always need to get a pointer to the first
 *  value in memory that should be written. The memory address of the
 *  z-coordinate of the first cell can be obtained using
 *  '&cells[0].m_coordinates[2]'. Therefore, in order to write your data you
 *  make the same call to writeArray() on both domains:
 *
 *      writeArray(&a_coordinates[2], "coordinates_2", 3);
 *
 *  What is the '3' doing there? Well, since the data values *in memory* are
 *  not contiguous (we only want to write the z-coordinates, i.e. we only need
 *  every third value), you need to provide a stride. The stride always refers
 *  to the layout in memory - in the file the values will always be contiguous.
 *
 *  Finally, when we check the file, we will see the following data in
 *  'coordinates_2':
 *
 *      | 0:z | 1:z | 2:z | 3:z | 4:z | ... | 998:z | 999:z |
 *
 *  Now, that wasn't so hard, was it?
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::writeArray(const T* array, const MString& name, size_type memoryStride,
                                         size_type diskStride) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  using namespace maia::parallel_io;

  // Prevent rewriting data in read-only mode
  if(m_fileMode == PIO_READ) {
    TERMM(1, "Cannot write data in read-only mode!");
  }

  // Abort if offsets were not set
  if(!m_offsetsAreSet) {
    TERMM(1, "Offsets have to be set before calling writeArray!");
  }

#ifdef MAIA_EXTRA_DEBUG
  // In debug mode, check if array was already written and give a warning
  if(m_writtenArrays.count(name) != 0) {
    std::cerr << "Warning: in " << AT_ << " (" << __FILE__ << ":" << __LINE__ << ") "
              << "the array '" << name << "' is written more than once. Make sure that this is the"
              << "intended behavior." << std::endl;
  }
  m_writtenArrays.insert(name);
  if(m_unwrittenArrays.count(name)) {
    m_unwrittenArrays.erase(name);
  }
#endif

  // Ensure reasonable parameter values
  if(memoryStride != -1 && memoryStride <= 0) {
    TERMM(1, "memoryStride must be greater than zero (or -1 to set automatically)!");
  }
  if(diskStride != -1 && diskStride <= 0) {
    TERMM(1, "diskStride must be greater than zero (or -1 to set automatically)!");
  }
  if(array == 0 && m_localCount[0] > 0) {
    TERMM(1, "Data pointer must point to a valid location (is: null pointer)!");
  }

  // Set actual stride if no explicit stride value was given
  size_type actualMemoryStride;
  if(memoryStride == -1) {
    actualMemoryStride = 1;
  } else {
    actualMemoryStride = memoryStride;
  }
  size_type actualDiskStride = (diskStride == -1) ? 1 : diskStride;

  // Save header (if not yet
  if(!m_isHeaderSaved) {
    derived().b_saveHeader();
    m_isHeaderSaved = true;
  }

  // Set number of chunks
  const size_type noChunks = (m_noChunks == -1) ? 1 : m_noChunks;

  size_type noDims = getDatasetNoDims(name);

  ScratchSpace<MPI_Offset> localCount(noDims, FUN_, "localCount");
  ScratchSpace<MPI_Offset> offset(noDims, FUN_, "offset");

  // Set counts for all dimensions except first to total dimension length
  // localCount[0] is set via setOffset(...)
  localCount[0] = m_localCount[0];
  for(MInt dimId = 1; dimId < noDims; dimId++) {
    localCount[dimId] = getArraySize(name, dimId);
  }
  for(MInt dimId = 0; dimId < noDims; dimId++) {
    offset[dimId] = m_offset[dimId];
  }

  // Write array
  derived().b_writeArray(array, name, noDims, &offset[0], &localCount[0], actualMemoryStride, noChunks,
                         actualDiskStride);
}


/**
 * \brief Write array data to file (multi-array version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the arrays (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data of the first array in memory.
 * \param[in] names Names of the arrays in the file.
 * \param[in] memoryStride Offset between two data values in *memory* (or -1 to set
 *                   automatically).
 * \param[in] diskStride Offset between two data values on *disk* (or -1 to set
 *                   automatically).
 *
 * The size and file offset have to be set using setOffset() before this method
 * is called.  If no stride is given (or -1), the stride will be set to the
 * number of provided names. A call to this method with n elements in 'names'
 * is *equivalent* to calling writeArray() n times, each time with the data
 * pointer 'array' being advanced by 1.
 *
 * Example:
 *
 * Suppose you would like to write all coordinates of all cells to a file, and
 * store them as separate variables 'coordinates_0' ... 'coordinates_2' (in
 * 3D). Similarly to the example in writeArray(), you have two domains with 550
 * cells on domain0, and 450 cells on domain1.
 *
 * First you have to call defineArray() and setOffset(). You call both with the
 * same arguments as described in writeArray(). The only difference is that
 * instead of defining each array individually, you can call
 * defineArray(maiabd_type, const vector<MString>&, size_type) with a vector
 * that holds your names (in our case the aforementioned 'coordinates_0',
 * 'coordinates_1', and 'coordinates_2'), which we call 'names' in this
 * example:
 *
 *     defineArray(PIO_FLOAT, names, 1000);
 *
 * This is *fully equivalent* to writing
 *
 *     defineArray(PIO_FLOAT, 'coordinates_0', 1000);
 *     defineArray(PIO_FLOAT, 'coordinates_1', 1000);
 *     defineArray(PIO_FLOAT, 'coordinates_2', 1000);
 *
 * Next, to set the offset, you call setOffset():
 *
 *     setOffset(550, 0);   // <-- This line on domain0. We have 550 cells and
 *                          //     start at the beginning.
 *     setOffset(450, 550); // <-- This line on domain1. We have 450 cells and
 *                          //     leave a gap of 550.
 *
 * Finally, we can write the data. Using the vector-of-strings version of
 * writeArray(maiabd_type, const vector<MString>&, size_type), we only need
 * one call:
 *
 *     writeArray(&a_coordinates[0], names, 3);
 *
 * This is *fully equivalent* to writing
 *
 *     writeArray(&a_coordinates[0], 'coordinates_0', 3);
 *     writeArray(&a_coordinates[1], 'coordinates_1', 3);
 *     writeArray(&a_coordinates[2], 'coordinates_2', 3);
 *
 * but more convenient :-). It gets even better: If you have want to write all
 * components of a multi-component array from memory to the file (as in this
 * example), writeArray(maiabd_type, const vector<MString>&, size_type) can do
 * the stride calculation automatically for you. So we could just leave it out
 * and write
 *
 *     writeArray(&cells[0].m_coordinates[0], names);
 *
 * instead.
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::writeArray(const T* array, const std::vector<MString>& names, size_type memoryStride,
                                         size_type diskStride) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  // Ensure reasonable parameter values
  if(memoryStride != -1 && memoryStride < static_cast<size_type>(names.size())) {
    TERMM(1, "memoryStride must be greater than or equal to number of "
             "variables (or -1 to set automatically)!");
  }
  if(diskStride != -1 && diskStride < static_cast<size_type>(names.size())) {
    TERMM(1, "diskStride must be greater than or equal to number of "
             "variables (or -1 to set automatically)!");
  }

  // Set actual stride if no explicit stride value was given
  size_type actualMemoryStride;
  if(memoryStride == -1) {
    actualMemoryStride = static_cast<size_type>(names.size());
  } else {
    actualMemoryStride = memoryStride;
  }
  size_type actualDiskStride = (diskStride == -1) ? 1 : diskStride;

  // Write arrays from vector using normal writeArray call
  for(std::vector<MString>::size_type i = 0; i < names.size(); i++) {
    writeArray(array + i, names[i], actualMemoryStride, actualDiskStride);
  }
}

/**
 * \brief Write multidimensional array data to file. <b>[MPI]</b>
 *
 * \author Marian Albers <m.albers@aia.rwth-aachen.de>
 * \date 2022-09-029
 *
 * \tparam T Data type of the array (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data in memory.
 * \param[in] path Group name in which the array shall be located (works only with HDF5).
 * \param[in] name Name of the array in the file.
 * \param[in] noDims Number of dimensions of the array
 * \param[in] offset Offset of the of the array in the dataset in each dimension, must be of size noDims
 * \param[in] localCount Size of the array in each dimension, must of of size noDims
 *
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::writeArray(const T* array, const MString& path, const MString& name,
                                         const size_type noDims, size_type* offset, size_type* localCount) {
  TRACE();

  ScratchSpace<size_type> ghost(noDims, FUN_, "ghost");
  ghost.fill(0);

  derived().b_writeArray(array, path, name, noDims, offset, localCount, &ghost[0]);
}


/**
 * \brief Write multidimensional array data to file [ghost cell version]. <b>[MPI]</b>
 *
 * \author Marian Albers <m.albers@aia.rwth-aachen.de>
 * \date 2022-09-029
 *
 * \tparam T Data type of the array (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data in memory.
 * \param[in] path Group name in which the array shall be located (works only with HDF5).
 * \param[in] name Name of the array in the file.
 * \param[in] noDims Number of dimensions of the array
 * \param[in] offset Offset of the of the array in the dataset in each dimension, must be of size noDims
 * \param[in] localCount Size of the array in each dimension, must of of size noDims
 * \param[in] ghost Number of ghost cells in memory preceding the actual values
 *
 * This method is specifically for storing two- or three-dimensional arrays
 * from the structured solver. Typically the in-memory array has a number of
 * ghost cells at the beginning in every dimension (and at the end) which
 * should be added to the offset when reading from memory, but should be ignored
 * for the dataset in the file that you're writing to. That is, the ghost values
 * are used for the access of the data in memory but ignored when they are being
 * written to the actual dataset in the file.
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::writeArray(const T* array,
                                         const MString& path,
                                         const MString& name,
                                         const size_type noDims,
                                         size_type* offset,
                                         size_type* localCount,
                                         size_type* ghost) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  using namespace maia::parallel_io;

  // Prevent rewriting data in read-only mode
  if(m_fileMode == PIO_READ) {
    TERMM(1, "Cannot write data in read-only mode!");
  }

#ifdef MAIA_EXTRA_DEBUG
  // In debug mode, check if array was already written and give a warning
  if(m_writtenArrays.count(name) != 0) {
    std::cerr << "Warning: in " << AT_ << " (" << __FILE__ << ":" << __LINE__ << ") "
              << "the array '" << name << "' is written more than once. Make sure that this is the"
              << "intended behavior." << std::endl;
  }
  m_writtenArrays.insert(name);
  if(m_unwrittenArrays.count(name)) {
    m_unwrittenArrays.erase(name);
  }
#endif

  // Save header (if not yet
  if(!m_isHeaderSaved) {
    derived().b_saveHeader();
    m_isHeaderSaved = true;
  }

  // Ensure reasonable parameter values
  if(array == 0 && localCount[0] > 0) {
    TERMM(1, "Data pointer must point to a valid location (is: null pointer)!");
  }

  // Write array
  derived().b_writeArray(array, path, name, noDims, offset, localCount, ghost);
}


/**
 * \brief Write scalar data to file. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the scalar (can be either MFloat or MInt and must
 *           match the previous definition of the scalar).
 * \param[in] scalar Value that should be written to file.
 * \param[in] name Name of the scalar in the file.
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::writeScalar(T scalar, const MString& name) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  using namespace maia::parallel_io;

  // Prevent rewriting data in read-only mode
  if(m_fileMode == PIO_READ) {
    mTerm(1, AT_, "Cannot write data in read-only mode!");
  }

#ifndef NDEBUG
  // In debug mode check if scalar value is the same on all ranks
  const T localValue = scalar;
  ScratchSpace<T> globalValues(m_noDomains, AT_, "globalValues");
  // Gather scalar values on rank 0
  MPI_Gather(&localValue, 1, maia::type_traits<T>::mpiType(), &globalValues[0], 1, maia::type_traits<T>::mpiType(), 0,
             m_mpiComm, AT_, "localValue", "globalValues[0]");

  // Compare values
  if(m_domainId == 0) {
    for(MInt i = 0; i < m_noDomains; i++) {
      MBool equal = true;
      // Note: cannot use == for floating point values
      if(globalValues[i] < localValue) {
        equal = false;
      }
      if(globalValues[i] > localValue) {
        equal = false;
      }
      if(!equal) {
        TERMM(1, "Error: scalar variable has not the same value on all domains, rank 0: " + std::to_string(localValue)
                     + " != " + std::to_string(globalValues[i]) + " (domain " + std::to_string(i) + ")");
      }
    }
  }
#endif

#ifdef MAIA_EXTRA_DEBUG
  // In extra debug mode, check if scalar was already written and give a warning
  if(m_writtenScalars.count(name) != 0) {
    std::cerr << "Warning: in " << AT_ << " (" << __FILE__ << ":" << __LINE__ << ") "
              << "the scalar '" << name << "' is written more than once. Make sure that this is the"
              << "intended behavior." << std::endl;
  }
  m_writtenScalars.insert(name);
  if(m_unwrittenScalars.count(name)) {
    m_unwrittenScalars.erase(name);
  }
#endif

  // Save header (if not yet done)
  if(!m_isHeaderSaved) {
    derived().b_saveHeader();
    m_isHeaderSaved = true;
  }

  derived().b_writeScalar(scalar, name);
}


/**
 * \brief Set a file or dataset attribute. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the attribute (can be either MFloat, MInt, or
 *           MString).
 * \param[in] value Value of the attribute.
 * \param[in] name Name of the attribute (must not be empty).
 * \param[in] datasetName Attaches the attribute to this dataset (i.e. array or
 *                        scalar). If empty, use as file attribute.
 *
 * You can add attributes to datasets as soon as they were defined with
 * defineArray()/ defineScalar(). You can also use this to overwrite an
 * attribute of the same name.  When you call this method before the first read
 * or write operation, there is no MPI overhead.
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::setAttribute(const T& value, const MString& name, const MString& datasetName) {
  TRACE();
  setAttributes(&value, name, 1, datasetName);
}

// Provide simple overload in order to accept string literals
template <class Backend>
void ParallelIoBase<Backend>::setAttribute(const MChar* value, const MString& name, const MString& datasetName) {
  const MString tmpValue = MString(value);
  setAttributes(&tmpValue, name, 1, datasetName);
}

template <class Backend>
template <class T>
void ParallelIoBase<Backend>::setAttributes(const T* value,
                                            const MString& name,
                                            size_type totalCount,
                                            const MString& datasetName) {
  TRACE();

#ifdef DISABLE_OUTPUT
  return;
#endif

  using namespace maia::parallel_io;

  // Prevent rewriting data in read-only mode
  if(m_fileMode == PIO_READ) {
    mTerm(1, AT_, "Cannot write data in read-only mode!");
  }

  // Ensure reasonable parameter values
  if(name.empty()) {
    mTerm(1, AT_, "Invalid name! Name must not be empty.");
  }
  ASSERT(totalCount > 0, "");

  derived().b_setAttribute(value, name, datasetName, totalCount);
}


//------------------------------------------------------------------------------
// Methods to read data or attributes
//------------------------------------------------------------------------------


/**
 * \brief Read array data from file. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the array (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data in memory.
 * \param[in] name Name of the array in the file.
 * \param[in] memoryStride Offset between two data values in *memory* (or -1 to set
 *                   automatically).
 * \param[in] diskStride Offset between two data values on *disk* (or -1 to set
 *                   automatically).
 *
 * The size and file offset have to be set using setOffset() before this method
 * is called.  If no stride is given (or -1), the stride will be set to 1.
 *
 * *Note:* This method is the exact counterpart to writeArray().
 *
 * Example:
 *
 * Suppose you would like to read in the z-coordinates you previously saved to
 * a file in the example of writeArray(maiabd_type, const MString&,
 * size_type). In this case, you just set the offset as before
 *
 *     setOffset(550, 0);   // <-- This line on domain0. We have 550 cells and
 *                          //     start at the beginning.
 *     setOffset(450, 550); // <-- This line on domain1. We have 450 cells and
 *                          //     leave a gap of 550.
 *
 * and then call readArray() with the exact same arguments as you used for
 * writeArray():
 *
 *     readArray(&a_coordinates[2], "coordinates_2", 3);
 *
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::readArray(T* array, const MString& name, size_type memoryStride, size_type diskStride) {
  TRACE();

  // Abort if offsets were not set
  if(!m_offsetsAreSet) {
    mTerm(1, AT_, "Offsets have to be set before calling readArray!");
  }

#ifdef MAIA_EXTRA_DEBUG
  // In debug mode, check if array was already read and give a warning
  if(m_readArrays.count(name) != 0) {
    std::cerr << "Warning: in " << AT_ << " (" << __FILE__ << ":" << __LINE__ << ") "
              << "the array '" << name << "' is read more than once. Make sure that this is the"
              << "intended behavior." << std::endl;
  }
  m_readArrays.insert(name);
#endif

  // Ensure reasonable parameter values
  if(memoryStride != -1 && memoryStride <= 0) {
    mTerm(1, AT_, "memoryStride must be greater than zero (or -1 to set automatically)!");
  }
  if(diskStride != -1 && diskStride <= 0) {
    mTerm(1, AT_, "diskStride must be greater than zero (or -1 to set automatically)!");
  }
  if(array == 0 && m_localCount[0] > 0) {
    mTerm(1, AT_, "Data pointer must point to a valid location (is: null pointer)!");
  }

  // Set actual stride if no explicit stride value was given
  size_type actualMemoryStride;
  if(memoryStride == -1) {
    actualMemoryStride = 1;
  } else {
    actualMemoryStride = memoryStride;
  }
  size_type actualDiskStride = (diskStride == -1) ? 1 : diskStride;

  // Set number of chunks
  const size_type noChunks = (m_noChunks == -1) ? 1 : m_noChunks;

  size_type noDims = getDatasetNoDims(name);

  ScratchSpace<MPI_Offset> localCount(noDims, FUN_, "localCount");
  ScratchSpace<MPI_Offset> offset(noDims, FUN_, "offset");

  localCount[0] = m_localCount[0];
  // Set counts for all dimensions except first to total dimension length
  // m_localCount[0] is set via setOffset(...)
  for(MInt dimId = 1; dimId < noDims; dimId++) {
    localCount[dimId] = getArraySize(name, dimId);
  }
  for(MInt dimId = 0; dimId < noDims; dimId++) {
    offset[dimId] = m_offset[dimId];
  }

  // Read array
  derived().b_readArray(array, name, noDims, &offset[0], &localCount[0], actualMemoryStride, noChunks,
                        actualDiskStride);
}

/**
 * \brief Read array data from file. <b>[MPI]</b>
 *
 * \author Marian Albers <m.albers@aia.rwth-aachen.de>
 * \date 2022-09-029
 *
 * \tparam T Data type of the array (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data in memory.
 * \param[in] name Name of the array in the file.
 * \param[in] path Group name in which the array shall be located (works only with HDF5).
 * \param[in] noDims Number of dimensions of the array
 * \param[in] offset Offset of the of the array in the dataset in each dimension, must be of size noDims
 * \param[in] localCount Size of the array in each dimension, must of of size noDims
 *
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::readArray(T* array,
                                        const MString& datasetName,
                                        const MString& name,
                                        const size_type noDims,
                                        const size_type* offset,
                                        const size_type* localCount) {
  TRACE();

#ifdef MAIA_EXTRA_DEBUG
  // In debug mode, check if array was already read and give a warning
  if(m_readArrays.count(name) != 0) {
    std::cerr << "Warning: in " << AT_ << " (" << __FILE__ << ":" << __LINE__ << ") "
              << "the array '" << name << "' is read more than once. Make sure that this is the"
              << "intended behavior." << std::endl;
  }
  m_readArrays.insert(name);
#endif

  if(array == 0 && localCount[0] > 0) {
    mTerm(1, AT_, "Data pointer must point to a valid location (is: null pointer)!");
  }

  // Read array
  derived().b_readArray(array, datasetName, name, noDims, offset, localCount);
}

/**
 * \brief Read array data from file [no offset version]. <b>[MPI]</b>
 *
 * \author Marian Albers <m.albers@aia.rwth-aachen.de>
 * \date 2022-09-029
 *
 * \tparam T Data type of the array (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data in memory.
 * \param[in] name Name of the array in the file.
 * \param[in] path Group name in which the array shall be located (works only with HDF5).
 * \param[in] noDims Number of dimensions of the array
 * \param[in] localCount Size of the array in each dimension, must of of size noDims
 *
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::readArray(T* array, const MString& datasetName, const MString& name,
                                        const size_type noDims, const size_type* localCount) {
  TRACE();

#ifdef MAIA_EXTRA_DEBUG
  // In debug mode, check if array was already read and give a warning
  if(m_readArrays.count(name) != 0) {
    std::cerr << "Warning: in " << AT_ << " (" << __FILE__ << ":" << __LINE__ << ") "
              << "the array '" << name << "' is read more than once. Make sure that this is the"
              << "intended behavior." << std::endl;
  }
  m_readArrays.insert(name);
#endif

  if(array == 0 && localCount[0] > 0) {
    mTerm(1, AT_, "Data pointer must point to a valid location (is: null pointer)!");
  }

  ScratchSpace<size_type> offset(noDims, FUN_, "offset");
  offset.fill(0);

  // Read array
  derived().b_readArray(array, datasetName, name, noDims, &offset[0], localCount);
}

/**
 * \brief Read array data from file (multi-array version). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the arrays (can be either MFloat or MInt and must
 *           match the previous definition of the array).
 * \param[in] array Pointer to the data of the first array in memory.
 * \param[in] names Names of the arrays in the file.
 * \param[in] memoryStride Offset between two data values in *memory* (or -1 to set
 *                   automatically).
 * \param[in] diskStride Offset between two data values on *disk* (or -1 to set
 *                   automatically).
 *
 *  The size and file offset have to be set using setOffset() before this
 *  method is called.  If no stride is given (or -1), the stride will be set to
 *  the number of provided names. A call to this method with n elements in
 *  'names' is *equivalent* to calling readArray() n times, each time with the
 *  data pointer 'array' being advanced by 1.
 *
 *  Example:
 *
 *  Suppose you would like to read in the all coordinates you previously saved
 *  to a file in the example of writeArray(maiabd_type, const
 *  vector<MString>&, size_type). In this case, you just set the offset as
 *  before
 *
 *      setOffset(550, 0);   // <-- This line on domain0. We have 550 cells and
 *                           //     start at the beginning.
 *      setOffset(450, 550); // <-- This line on domain1. We have 450 cells and
 *                           //     leave a gap of 550.
 *
 *  and then call readArray(maiabd_type, const vector<MString>&, size_type)
 *  with the exact same arguments as you used for writeArray(maiabd_type, const
 *  vector<MString>&, size_type):
 *
 *      readArray(&cells[0].m_coordinates[0], names);
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::readArray(T* array, const std::vector<MString>& names, size_type memoryStride,
                                        size_type diskStride) {
  TRACE();

  // Ensure reasonable parameter values
  if(memoryStride != -1 && memoryStride < static_cast<size_type>(names.size())) {
    mTerm(1, AT_,
          "memoryStride must be greater than or equal to number of "
          "variables (or -1 to set automatically)!");
  }
  if(diskStride != -1 && diskStride < static_cast<size_type>(names.size())) {
    mTerm(1, AT_,
          "diskStride must be greater than or equal to number of "
          "variables (or -1 to set automatically)!");
  }

  // Set actual stride if no explicit stride value was given
  size_type actualMemoryStride;
  if(memoryStride == -1) {
    actualMemoryStride = static_cast<size_type>(names.size());
  } else {
    actualMemoryStride = memoryStride;
  }
  size_type actualDiskStride = (diskStride == -1) ? 1 : diskStride;

  // Read arrays from vector using normal readArray call
  for(std::vector<MString>::size_type i = 0; i < names.size(); i++) {
    readArray(array + i, names[i], actualMemoryStride, actualDiskStride);
  }
}


/**
 * \brief Read scalar data from file. <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the scalar (can be either MFloat or MInt and must
 *           match the previous definition of the scalar).
 * \param[in] scalar Value that should be read from file.
 * \param[in] name Name of the scalar in the file.
 *
 * \details *Note:* This method is the exact counterpart to writeScalar().
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::readScalar(T* scalar, const MString& name) {
#ifdef MAIA_EXTRA_DEBUG
  // In debug mode, check if scalar was already read and give a warning
  if(m_readScalars.count(name) != 0) {
    std::cerr << "Warning: in " << AT_ << " (" << __FILE__ << ":" << __LINE__ << ") "
              << "the scalar '" << name << "' is read more than once. Make sure that this is the"
              << "intended behavior." << std::endl;
  }
  m_readScalars.insert(name);
#endif

  derived().b_readScalar(scalar, name);
}


/**
 * \brief Retrieve a file or dataset attribute.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \tparam T Data type of the attribute (can be either MFloat, MInt, or
 *           MString).
 * \param[in] value Value of the attribute.
 * \param[in] name Name of the attribute.
 * \param[in] datasetName Reads the attribute from this dataset (i.e. array or
 *                        scalar). If empty, retrieve as file attribute.
 */
template <class Backend>
template <class T>
void ParallelIoBase<Backend>::getAttribute(T* value, const MString& name, const MString& datasetName) {
  TRACE();
  getAttribute(value, name, 1, datasetName);
}

template <class Backend>
template <class T>
void ParallelIoBase<Backend>::getAttribute(T* value, const MString& name, const size_type totalCount,
                                           const MString& datasetName) {
  TRACE();
  derived().b_getAttribute(value, name, datasetName, totalCount);
}

template <class Backend>
void ParallelIoBase<Backend>::getAttributeCount(const MString& name, size_type* totalCount,
                                                const MString& datasetName) {
  TRACE();
  derived().b_getAttributeCount(name, totalCount, datasetName);
}

//------------------------------------------------------------------------------
// Methods to calculate and set offset
//------------------------------------------------------------------------------


/**
 * \brief Set the local and global counts, as well the local offset for array
 *        operations.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] localCount Number of elements in an array on the current domain.
 * \param[in] offset Starting offset in the *file*, from which this domain
 *                   should start writing/reading.
 * \param[in] noDims number of dimensions of the array
 *            for more than one dimension localCount and offset refer to the
 *            first dimension, for all other dimensions the count will be set
 *            to the total dimension length, multidimensional counts and
 *            offsets are not yet supported
 * \param[in] noChunks Number of Chunks for chunked IO.
 *
 * You have to call this method at least once before your first call to
 * readArray()/writeArray(). The reason for storing the offset and count
 * information with the object is that there are usually only a very limited
 * number of different array sizes in each file, and thus it is more convenient
 * to set the offsets once and write/read multiple times afterwards.
 *
 * The offset may be calculated by hand (or known beforehand), or may be
 * conveniently calculated using calcOffset().
 */
template <class Backend>
void ParallelIoBase<Backend>::setOffset(const size_type localCount,
                                        const size_type offset,
                                        const size_type noDims,
                                        const size_type noChunks) {
  TRACE();

  // Ensure reasonable parameter values
  if(localCount < 0) {
    TERMM(1, "The local size may not be less than zero (was:" + std::to_string(localCount) + ")!");
  }
  if(offset < 0) {
    TERMM(1, "The offset may not be less than zero (was:" + std::to_string(offset) + ")!");
  }
  if(noDims <= 0) {
    TERMM(1, "The number of dimensions must at least be one!");
  }
  if(noChunks == 0) {
    TERMM(1, "The number of chunks may not be zero!");
  }

  // TODO labels:IO add support for setting a multi-dimensional count/offset for reading/writing data

  // Set member variables
  m_localCount.reserve(noDims);
  m_localCount.assign(noDims, -1); // -1 as a precaution
  m_localCount[0] = localCount;

  // Note:
  // if noDims>1: {m_localCount[i], i>0} will be set to the total count for that
  // dimension before reading/writing the multi-d array
  // setting the local count for multiple dimensions is not yet supported
  m_offset.reserve(noDims);
  m_offset.assign(noDims, 0);
  m_offset[0] = offset;

  if(noChunks > 0) {
    m_noChunks = noChunks;
  }

  // Set state variable
  m_offsetsAreSet = true;
}


/**
 * \brief Set the local and global counts, as well the local offset for array
 *        operations (overload).
 */
template <class Backend>
void ParallelIoBase<Backend>::setOffset(const size_type localCount, const size_type offset, const size_type noDims) {
  TRACE();

  setOffset(localCount, offset, noDims, -1);
}


/**
 * \brief Set the local and global counts, as well the local offset for array
 *        operations (overload).
 */
template <class Backend>
void ParallelIoBase<Backend>::setOffset(const size_type localCount, const size_type offset) {
  TRACE();

  setOffset(localCount, offset, 1, -1);
}


/**
 * \brief Calculates the totalCount and offset information needed in
 *        setOffset(). <b>[MPI]</b>
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Konstantin Froehlich
 * \date 2013-04-09
 *
 * \param[in] localCount Number of elements in an array on the current domain.
 * \param[out] offset Starting offset in the *file*, from which this domain
 *                    should start writing/reading (is zero if localCount=0).
 * \param[out] totalCount Total number of elements in an array across all
 *                        domains.
 * \param[in] mpiComm MPI communicator for which the totalCount and the offset
 *                    should be calculated.
 * \param[in] noChunks Number of Chunks for chunked IO.
 *
 * \details This must be a collective call from all ranks in mpiComm.
 */
template <class Backend>
void ParallelIoBase<Backend>::calcOffset(size_type localCount, size_type* offset, size_type* totalCount,
                                         const MPI_Comm& mpiComm, size_type* noChunks) {
  TRACE();

  // Ensure reasonable parameter values
  if(localCount < 0) {
    TERMM(1, "The local size may not be less than zero (was:" + std::to_string(localCount) + ")!");
  }

  // Get rank and size
  MInt domainId, noDomains;
  MPI_Comm_rank(mpiComm, &domainId);
  MPI_Comm_size(mpiComm, &noDomains);

  // Calculate offset and total size values
  MPI_Exscan(&localCount, offset, 1, s_mpiSizeType, MPI_SUM, mpiComm, AT_, "localCount", "offset");
  if(domainId == 0) {
    offset[0] = 0;
  }

  totalCount[0] = offset[0] + localCount;
  MPI_Bcast(totalCount, 1, s_mpiSizeType, (noDomains - 1), mpiComm, AT_, "totalCount");

  // If a domain has zero size set the offset to zero such that reading/writing
  // can not result in an error caused by exceeding the array size
  if(localCount == 0) {
    offset[0] = 0;
  }

  // If noChunks is not a null pointer, calculate & set number of chunks
  if(noChunks != nullptr) {
    // Divide total count by maximum chunk size
    noChunks[0] = totalCount[0] / s_maxChunkSize;

    // If there is a remainder for the division, increment number by one
    if(totalCount[0] % s_maxChunkSize > 0) {
      noChunks[0] += 1;
    }
  }
}


/// \brief Auxiliary method that aborts MAIA if type is not a valid type constant for ParallelIo.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2017-06-23
///
/// \param[in] type The type constant (integer) to check for validity.
template <class Backend>
void ParallelIoBase<Backend>::validateType(const maiabd_type type) const {
  TRACE();

  using namespace maia::parallel_io;

  if(type != PIO_FLOAT && type != PIO_INT && type != PIO_LONG && type != PIO_STRING && type != PIO_UCHAR
     && type != PIO_ULONGLONG) {
    TERMM(1, "Invalid data type! Must be PIO_FLOAT or PIO_INT or PIO_LONG or PIO_STRING or PIO_UCHAR.");
  }
}


#endif /* PARALLELIO_H */
