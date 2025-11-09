// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// NOTE: commented on purpose to prevent multiple includes of mpi.h as it should only be included once in mpioverride.h
// which is the header that should be used everywhere
/* #ifndef MAIA_MPI_H_ */
/* #define MAIA_MPI_H_ */
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file wraps the mpi implementation
///
/// Note: this allows to disable warnings that are outside our control
////////////////////////////////////////////////////////////////////////////////
#include "INCLUDE/maiatypes.h"
#include "compiler_config.h"
#include "config.h"
////////////////////////////////////////////////////////////////////////////////
/// \brief ignore deprecated c++ binding of mpi
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
/// \brief avoid conflict between mpi in c++ binding and stdio.h
#ifndef MPICH_IGNORE_CXX_SEEK
#define MPICH_IGNORE_CXX_SEEK
#endif
////////////////////////////////////////////////////////////////////////////////
/// \brief Include the mpi.h header
#if defined(MAIA_GCC_COMPILER)
/// \brief suppress pedantic and long-long warnings if using MAIA_GCC_COMPILER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wlong-long"
/// Suppress more warnings on Hawk
#if defined(HOST_Hawk)
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wredundant-decls"
#endif
#include <mpi.h>
#if defined(HOST_Hawk)
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#endif
#pragma GCC diagnostic pop
#elif defined(MAIA_CLANG_COMPILER)
/// Suppress cast alignment warnings if using clang
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#include <mpi.h>
#pragma clang diagnostic pop
#else
#include <mpi.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#if MPI_VERSION < 3 && !defined(HOST_Klogin)
#define MPI_Iallreduce MPIX_Iallreduce
#endif

// minimum upper bound for MPI tags as defined by the MPI standard
#ifndef MPI_TAG_UB
#define MPI_TAG_UB 32767
#endif

/// Define an MPI null request with MPI_Request datatype to avoid templating errors in mAlloc when
/// MPI_REQUEST_NULL is passed as a default value
static const MPI_Request MPI_REQ_NULL = MPI_REQUEST_NULL;

/* #endif */
