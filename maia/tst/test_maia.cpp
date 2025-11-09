// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

// This is only needed, as mEnvironment is declared 'extern' in
// src/UTIL/functions.cpp
#include "environment.h"
Environment* mEnvironment = nullptr;
