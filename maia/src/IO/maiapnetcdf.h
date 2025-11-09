// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2019 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier:    LGPL-3.0-only


#ifndef MAIA_PNETCDF_H_
#define MAIA_PNETCDF_H_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file wraps the pnetcdf implementation
///
/// Note: this allows to disable warnings that are outside our control
////////////////////////////////////////////////////////////////////////////////
#include "compiler_config.h"
#include "config.h"
////////////////////////////////////////////////////////////////////////////////
/// \brief Include the pnetcdf.h header
#ifdef MAIA_GCC_COMPILER
/// \brief suppress long-long warning if using MAIA_GCC_COMPILER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wlong-long"
#include <pnetcdf.h>
#pragma GCC diagnostic pop
#elif not defined(MAIA_MS_COMPILER)
#include <pnetcdf.h>
#endif

#endif
