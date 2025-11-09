// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PARALLELIO_HDF5_H
#define PARALLELIO_HDF5_H

#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic ignored "-Wredundant-tags"
#endif
#include "hdf5.h"
#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic pop
#endif
#include "parallelio.h"

class ParallelIoHdf5 : public ParallelIoBase<ParallelIoHdf5> {
  friend ParallelIoBase<ParallelIoHdf5>;

  // Public member functions
 public:
  // Constructor
  ParallelIoHdf5(const MString& fileName, MInt fileMode, const MPI_Comm& mpiComm);
  // Destructor
  ~ParallelIoHdf5() override;

  // Hdf5-specific methods (private)
 private:
  // Static file system-related methods
  static MBool b_isValidFile(const MString& name, const MPI_Comm& mpiComm);
  static const MString b_fileExt();

  // File methods
  void b_open();
  void b_saveHeader();
  void b_addAdditionalHeader();
  void b_writeAdditionalData();


  // Define mode methods
  void b_defineArray(maiabd_type type, const MString& name, size_type noDims, size_type* totalCount);
  void b_defineArray(maiabd_type, const MString&, const MString&, const size_type, const size_type*);
  void b_defineScalar(maiabd_type type, const MString& name);

  void b_createGroup(const MString& path);
  size_type b_getNoDatasets(const MString& path);
  size_type b_getNoGroups(const MString& path);

  // Inquiry methods
  MBool b_hasDataset(const MString& name, const size_type dimension);
  MBool b_hasDataset(const MString& name, const MString& path);
  MBool b_hasObject(const MString& path);
  MInt b_getDatasetType(const MString& name);
  void b_getDatasetNames(std::vector<MString>& names, const size_type dimension);
  void b_getDatasetNames(std::vector<MString>& names, const MString& path);
  void b_getGroupNames(std::vector<MString>& names, const MString& path);
  size_type b_getDatasetNoDims(const MString& name);
  size_type b_getDatasetNoDims(const MString& path, const MString& name);
  size_type b_getDatasetSize(const MString& name, const size_type dimensionId = 0);
  void b_getDatasetSize(const MString& name, const MString& path, size_type noDims, size_type* data);
  MBool b_hasAttribute(const MString& name, const MString& path = "");
  MInt b_getAttributeType(const MString& name, const MString& datasetName = "");

  // Data mode methods
  template <class T>
  void b_writeArray(const T* array, const MString& name, const size_type noDims, MPI_Offset* start, MPI_Offset* count,
                    MPI_Offset memoryStride, const size_type noChunks, MPI_Offset diskStride);
  template <class T>
  void b_writeArray(const T*, const MString&, const MString&, const size_type, const size_type*, const size_type*,
                    const size_type*);
  template <class T>
  void b_writeScalar(T scalar, const MString& name);

  template <class T>
  void b_readArray(T*, const MString&, const size_type, MPI_Offset*, MPI_Offset*, MPI_Offset, const size_type,
                   MPI_Offset);
  template <class T>
  void b_readArray(T*, const MString&, const MString&, const size_type, const size_type*, const size_type*);

  template <class T>
  void b_readScalar(T* scalar, const MString& name);

  // Attribute methods
  template <class T>
  void b_setAttribute(const T* value, const MString& name, const MString& datasetName = "",
                      const size_type totalCount = 1);
  template <class T>
  void b_createAttribute(const T* value, const MString& name, const MString& datasetName, hid_t type_id,
                         const size_type totalCount);
  template <class T>
  void b_getAttribute(T* value, const MString& name, const MString& datasetName = "", const size_type totalCount = 1);
  void b_getAttributeCount(const MString& name, size_type* totalCount, const MString& datasetName = "");

  // Auxiliary methods
  // b_error is static because it is used in b_isValidFile (also static)
  static void b_error(MInt status, const MString& name, const MString& location);
  static void b_warning(MInt status, const MString& name, const MString& location);

  // Hdf5-related member variables
 private:
  hid_t b_h5Id = -1;
  // member to manage the I/O operations (collectively)
  hid_t b_h5DatasetXferHandle;
  hid_t b_h5FileXferHandle;
};
#endif /* PARALLELIO_H */
