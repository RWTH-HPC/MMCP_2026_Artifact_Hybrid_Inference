// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "couplingdgape.h"

const std::array<MString, MeanVariables<2>::totalNoMeanVars> MeanVariables<2>::names = {
    {"wxv_x",     "wxv_y",     "um",          "vm",          "rhom",          "pm",
     "c0",        "dc0_dx",    "dc0_dy",      "vort_z",      "gradm_u_x",     "gradm_u_y",
     "gradm_v_x", "gradm_v_y", "gradm_rho_x", "gradm_rho_y", "gradm_p_x",     "gradm_p_y",
     "rhodivu_x", "rhodivu_y", "ugradrho_x",  "ugradrho_y",  "mgrad_p_rho_x", "mgrad_p_rho_y",
     "ugradu_x",  "ugradu_y"}};
const std::array<MString, MeanVariables<3>::totalNoMeanVars> MeanVariables<3>::names = {{"wxv_x",
                                                                                         "wxv_y",
                                                                                         "wxv_z",
                                                                                         "um",
                                                                                         "vm",
                                                                                         "wm",
                                                                                         "rhom",
                                                                                         "pm",
                                                                                         "c0",
                                                                                         "dc0_dx",
                                                                                         "dc0_dy",
                                                                                         "dc0_dz",
                                                                                         "vort_x",
                                                                                         "vort_y",
                                                                                         "vort_z",
                                                                                         "gradm_u_x",
                                                                                         "gradm_u_y",
                                                                                         "gradm_u_z",
                                                                                         "gradm_v_x",
                                                                                         "gradm_v_y",
                                                                                         "gradm_v_z",
                                                                                         "gradm_w_x",
                                                                                         "gradm_w_y",
                                                                                         "gradm_w_z",
                                                                                         "gradm_rho_x",
                                                                                         "gradm_rho_y",
                                                                                         "gradm_rho_z",
                                                                                         "gradm_p_x",
                                                                                         "gradm_p_y",
                                                                                         "gradm_p_z",
                                                                                         "rhodivu_x",
                                                                                         "rhodivu_y",
                                                                                         "rhodivu_z",
                                                                                         "ugradrho_x",
                                                                                         "ugradrho_y",
                                                                                         "ugradrho_z",
                                                                                         "mgrad_p_rho_x",
                                                                                         "mgrad_p_rho_y",
                                                                                         "mgrad_p_rho_z",
                                                                                         "ugradu_x",
                                                                                         "ugradu_y",
                                                                                         "ugradu_z"}};

const std::array<MString, 6> MeanVariables<2>::nodeVarNames = {{"um", "vm", "rhom", "c0", "dc0_dx", "dc0_dy"}};
const std::array<MString, 8> MeanVariables<3>::nodeVarNames = {
    {"um", "vm", "wm", "rhom", "c0", "dc0_dx", "dc0_dy", "dc0_dz"}};
