// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "globalmpiinfo.h"
#include "mpioverride.h"

GlobalMpiInformation g_mpiInformation;


/// \brief Print all information of the given MPI_Info object.
void printMpiInfo(MPI_Info& mpiInfo) {
  MInt i, nkeys;

  MPI_Info_get_nkeys(mpiInfo, &nkeys, AT_);
  std::cerr << "MPI Info: nkeys = " << nkeys << std::endl;
  for(i = 0; i < nkeys; i++) {
    char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
    MInt valuelen, flag;

    MPI_Info_get_nthkey(mpiInfo, i, key, AT_);
    MPI_Info_get_valuelen(mpiInfo, key, &valuelen, &flag, AT_);
    MPI_Info_get(mpiInfo, key, valuelen + 1, value, &flag, AT_);
    std::cerr << "MPI Info: [" << i << "] key = " << key << ", value = " << value << std::endl;
  }
}
