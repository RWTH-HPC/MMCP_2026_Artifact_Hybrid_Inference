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


#ifndef MACRO_H_
#define MACRO_H_

#include <string>

////////////////////////////////////////////////////////////////////////////////
/// \file \brief All global macro definitions that are not compiler-specific
///              (and thus should go into compiler_config.h), that are not config-
///              related (and thus should go into config.h), and that are not
///              file/class-specific (e.g. debug macros should be in debug.h)
///              are defined here.
////////////////////////////////////////////////////////////////////////////////


/// Define macros to stringify literal and expanded macro arguments
///
/// STRINGIFY() can be used to stringify literal macro arguments (e.g.
/// STRINGIFY(__FILE__) becomes "__FILE__"), while XSTRINGIFY() will expand a
/// macro first (e.g. XSTRINGIFY(__FILE__) becomes "macro.h").
#define STRINGIFY(s) #s
#define XSTRINGIFY(s) STRINGIFY(s)

/// Define a short-hand macros for the location in the code (<file>:<line>)
#define LOC_ __FILE__ ":" XSTRINGIFY(__LINE__)


/// Define FUN_ for Paraview plugin
#ifdef PVPLUGIN
#ifndef FUN_
#define FUN_ __PRETTY_FUNCTION__
#endif
#endif

/// AT_ is used e.g. in the scratch and properties class to
/// provide informations about the calling function. Since the macro
/// __PRETTY_FUNCTION__ is no C++ standard, the macro FUN_ is used in MAIA, which
/// is set in compiler_config.h. Furthermore, the function name is augmented with
/// the location in the code (e.g. "void set(int) (maia.h:31)")
///
/// Note: __CALLING__FUNCTION__ is deprecated, since identifiers starting with
///       underscores are reserved for the compiler. Use FUN_ or AT_ instead.
#define AT_ std::string(FUN_) + " (" + LOC_ + ")"

/// Use this macro to indicate a non-used function parameter.
/// \details This is the preferred method for example in overloaded functions, where the
///          default implementation is empty. In that case the compiler would produce
///          a 'parameter not used' warning for each unused parameter. Using this macro
///          is better than removing the parameter name, as it allows the source code
///          reader to directly understand the meaning of each parameter and that it will
///          not be used in the current implementation.
#define NotUsed(a)

#endif // ifndef MACRO_H_
