// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "debug.h"

using namespace std;

MInt MDebug::m_debugLevel;
MBool MDebug::m_debugOn;
basic_string<char> MDebug::m_traceSpaces;
