// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_GLOBALMPIINFO_H
#define MAIA_GLOBALMPIINFO_H

#include <iostream>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "mpioverride.h"

/// Print all information of given MPI_Info object
void printMpiInfo(MPI_Info& mpiInfo);

////////////////////////////////////////////////////////////////////////////////
// Accessors and storage for global MPI information

/// Class to store global MPI information and to prevent accidental changes
class GlobalMpiInformation {
 public:
  void init(const MInt domainId, const MInt noDomains, const MPI_Comm maiaCommWorld) {
    m_globalDomainId = domainId;
    m_globalNoDomains = noDomains;
    m_maiaCommWorld = maiaCommWorld;

    initMPIInformation();
  }

  void printInfo() {
#ifdef MPI_IO_PRINT_INFO
    // Print MPI information on global rank 0
    if(m_globalDomainId == 0) {
      printMpiInfo(m_mpiInfo);
    }
#endif
  }

  MPI_Comm getMaiaCommWorld(){ return m_maiaCommWorld; }

 private:
  void initMPIInformation() {
    MPI_Info_create(&m_mpiInfo, AT_);

    // Set header align size to 10KB for netCDF files. Allows to append header data without the need
    // to move all variable data if the header size is exceeded (which may cause MPI I/O errors).
    // Source: https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/VariableAlignment
    MPI_Info_set(m_mpiInfo, "nc_header_align_size", "10240");
    // Note: possibility to set variable align size
    /* MPI_Info_set(m_mpiInfo, "nc_var_align_size", "4194304"); */

#if !defined(WITH_HDF5) && defined(MPI_IO_OPT) && defined(HOST_HAZELHEN)
    // taken from Cray Wiki: https://wickie.hlrs.de/platforms/index.php/MPI-IO,
    // see also PNetcdf documentation: http://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/HintsForPnetcdf
    if(m_globalNoDomains > 256) {
      MPI_Info_set(m_mpiInfo, (char*)"cb_align", (char*)"2");             /* Default: OMPI: none, CrayMPT: 2 */
      MPI_Info_set(m_mpiInfo, (char*)"cb_nodes_list", (char*)"*:*");      /* Default: OMPI: *:1, CrayMPT: *:* */
      MPI_Info_set(m_mpiInfo, (char*)"direct_io", (char*)"false");        /* Default: OMPI: none, CrayMPT: false */
      MPI_Info_set(m_mpiInfo, (char*)"romio_ds_read", (char*)"disable");  /* Default: OMPI: none, CrayMPT: disable */
      MPI_Info_set(m_mpiInfo, (char*)"romio_ds_write", (char*)"disable"); /* Default: OMPI: none, CrayMPT: disable */
      /* Let's reduce the number of aggregators, should be roughly 2 to 4 times the stripe-factor */
      // MPI_Info_set (m_mpiInfo, (char*)"cb_nodes", (char*)"8");
      /* Default: OMPI: set automatically to the number of distinct nodes; However TOO High */

      MPI_Info_set(m_mpiInfo, (char*)"ind_wr_buffer_size", (char*)"16777216");
      /* proposed by PNetcdf documentation */
      MPI_Info_set(m_mpiInfo, (char*)"striping_factor", (char*)"64");
      /* no. of I/O devices across which the file should be striped */
      MPI_Info_set(m_mpiInfo, (char*)"cb_nodes", (char*)"128");
    }
#endif

#if defined(MPI_IO_OPT) && defined(HOST_Hawk)
    if(m_globalNoDomains > 10000) { // TODO labels:HAWK,IO
      // NOTE: PNetcdf memory issue for large scale simulations. During the pnetcdf write call a significant amount of
      // memory is allocated (at least on Hawk; scales linear with noDomains), which is not freed thereafter.
      // Setting these romio hints solves the memory allocation problem, however it is not clear if this was responsible
      // for some incomplete written data files.
      // To be able to check for erroneous files you can enabled the fill mode for PNetcdf in config.h with
      // MAIA_NCMPI_FILL_VARIABLES = true
      // and check your files for fill values in the data (which should have been overwritten).
      MPI_Info_set(m_mpiInfo, (char*)"romio_cb_read", (char*)"disable");
      MPI_Info_set(m_mpiInfo, (char*)"romio_cb_write", (char*)"disable");
      if(m_globalDomainId == 0) {
        std::cerr << std::endl
                  << std::endl
                  << "NOTE: disabling ROMIO hints romio_cb_read/write to avoid PNetcdf/Hdf5 memory allocation issues "
                     "on HAWK... "
                  << std::endl
                  << "NOTE: see comment at " << AT_ << std::endl
                  << "NOTE: undefine MPI_IO_OPT to turn off the ROMIO hint changes." << std::endl
                  << std::endl
                  << std::endl;
      }

      // TODO labels:HAWK,IO check if it makes sense to disable these
      // MPI_Info_set(m_mpiInfo, (char*)"romio_ds_read", (char*)"disable");
      // MPI_Info_set(m_mpiInfo, (char*)"romio_ds_write", (char*)"disable");
    }
#endif

#ifdef MPI_IO_PRINT_INFO
    // Print MPI information on global rank 0
    if(m_globalDomainId == 0) {
      std::cerr << std::endl << "Global MPI information" << std::endl;
      printMpiInfo(m_mpiInfo);
    }
#endif
  }

  friend MInt globalDomainId();
  friend MInt globalNoDomains();
  friend const MPI_Info& globalMpiInfo();
  friend MPI_Comm& globalMaiaCommWorld();

  MInt m_globalDomainId = 0;
  MInt m_globalNoDomains = 1;
  MPI_Comm m_maiaCommWorld = MPI_COMM_NULL;
  MPI_Info m_mpiInfo = MPI_INFO_NULL;
};

extern GlobalMpiInformation g_mpiInformation;

/// Return global domain id
inline MInt globalDomainId() { return g_mpiInformation.m_globalDomainId; }
/// Return global number of domains
inline MInt globalNoDomains() { return g_mpiInformation.m_globalNoDomains; }
/// Return global MPI information
inline const MPI_Info& globalMpiInfo() { return g_mpiInformation.m_mpiInfo; }

inline MPI_Comm& globalMaiaCommWorld() { return g_mpiInformation.m_maiaCommWorld; }

#endif // MAIA_GLOBALMPIINFO_H
