// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "enums.h"

#include <iostream>

using namespace std;

/** \brief This global function translates strings in their corresponding
 *  enum values (integer values).
 *  (Maybe something like namespaces should be introduced to avoid
 *   name conflicts)
 *
 */

MInt string2enum(MString theString) {
  if(theString == MString("NETCDF")) return (MInt)NETCDF;
  if(theString == MString("HDF")) return (MInt)HDF;
  if(theString == MString("TOML")) return (MInt)TOML;
  if(theString == MString("ASCII")) return (MInt)ASCII;
  if(theString == MString("BINARY")) return (MInt)BINARY;
  if(theString == MString("VTK")) return (MInt)VTK;
  if(theString == MString("VTU")) return (MInt)VTU;
  if(theString == MString("MAIA_STRUCTURED")) return (MInt)MAIA_STRUCTURED;
  if(theString == MString("MAIA_FINITE_VOLUME")) return (MInt)MAIA_FINITE_VOLUME;
  if(theString == MString("MAIA_FV_APE")) return (MInt)MAIA_FV_APE;
  if(theString == MString("MAIA_FV_GEQU_PV")) return (MInt)MAIA_LS_COMBUSTION;
  if(theString == MString("MAIA_LEVELSET")) return (MInt)MAIA_LEVELSET;
  if(theString == MString("MAIA_LS_SOLVER")) return (MInt)MAIA_LS_SOLVER;
  if(theString == MString("MAIA_LEVELSET_SOLVER")) return (MInt)MAIA_LEVELSET_SOLVER;
  if(theString == MString("MAIA_POST_DATA")) return (MInt)MAIA_POST_DATA;
  if(theString == MString("MAIA_LS_FV")) return (MInt)MAIA_LS_FV;
  if(theString == MString("MAIA_LS_COMBUSTION")) return (MInt)MAIA_LS_COMBUSTION;
  if(theString == MString("MAIA_MULTI_LS")) return (MInt)MAIA_MULTI_LS;
  if(theString == MString("MAIA_FV_LEVELSET")) return (MInt)MAIA_FV_LEVELSET;
  if(theString == MString("MAIA_LATTICE_BOLTZMANN")) return (MInt)MAIA_LATTICE_BOLTZMANN;
  if(theString == MString("MAIA_FINITE_CELL")) return (MInt)MAIA_FINITE_CELL;
  if(theString == MString("MAIA_FV_MB")) return (MInt)MAIA_FV_MB;
  if(theString == MString("MAIA_FV_MB_NEW_RK")) return (MInt)MAIA_FV_MB_NEW_RK;
  if(theString == MString("MAIA_DISCONTINUOUS_GALERKIN")) return (MInt)MAIA_DISCONTINUOUS_GALERKIN;
  if(theString == MString("MAIA_DG_MULTISOLVER")) return (MInt)MAIA_DG_MULTISOLVER;
  if(theString == MString("MAIA_DG_FV")) return (MInt)MAIA_DG_FV;
  if(theString == MString("MAIA_RUNGE_KUTTA")) return (MInt)MAIA_RUNGE_KUTTA;
  if(theString == MString("MAIA_RUNGE_KUTTA_STRUCTURED")) return (MInt)MAIA_RUNGE_KUTTA_STRUCTURED;
  if(theString == MString("MAIA_RUNGE_KUTTA_GEQU_PV")) return (MInt)MAIA_RUNGE_KUTTA_GEQU_PV;
  if(theString == MString("MAIA_RUNGE_KUTTA_LEVELSET")) return (MInt)MAIA_RUNGE_KUTTA_LEVELSET;
  if(theString == MString("MAIA_SEMI_LAGRANGE_LEVELSET")) return (MInt)MAIA_SEMI_LAGRANGE_LEVELSET;
  if(theString == MString("MAIA_SEMI_LAGRANGE_LEVELSET_LB")) return (MInt)MAIA_SEMI_LAGRANGE_LEVELSET_LB;
  if(theString == MString("MAIA_ANALYTIC_LEVELSET")) return (int)MAIA_ANALYTIC_LEVELSET;
  if(theString == MString("MAIA_LEVELSET_SURFACE")) return (int)MAIA_LEVELSET_SURFACE;
  if(theString == MString("MAIA_RUNGE_KUTTA_MB_SEMI_LAGRANGE_LEVELSET"))
    return (MInt)MAIA_RUNGE_KUTTA_MB_SEMI_LAGRANGE_LEVELSET;
  if(theString == MString("MAIA_RUNGE_KUTTA_MB_LEVELSET")) return (MInt)MAIA_RUNGE_KUTTA_MB_LEVELSET;
  if(theString == MString("MAIA_UNIFIED")) return (MInt)MAIA_UNIFIED;
  if(theString == MString("MAIA_PARTICLE")) return (MInt)MAIA_PARTICLE;
  if(theString == MString("MAIA_ACOUSTIC_ANALOGY")) return (MInt)MAIA_ACOUSTIC_ANALOGY;
  if(theString == MString("MAIA_RIGID_BODIES")) return (MInt)MAIA_RIGID_BODIES;

  // Grid types
  if(theString == MString("MAIA_GRID_CARTESIAN")) return (MInt)MAIA_GRID_CARTESIAN;
  if(theString == MString("MAIA_GRID_STRUCTURED")) return (MInt)MAIA_GRID_STRUCTURED;
  if(theString == MString("MAIA_GRID_NONE")) return (MInt)MAIA_GRID_NONE;

  // sensors
  if(theString == MString("DERIVATIVE")) return (MInt)DERIVATIVE;
  if(theString == MString("ENTROPY_GRADIENT")) return (MInt)ENTROPY_GRADIENT;
  if(theString == MString("ENTROPY_QUOTIENT")) return (MInt)ENTROPY_QUOTIENT;
  if(theString == MString("VORTICITY")) return (MInt)VORTICITY;
  if(theString == MString("INTERFACE")) return (MInt)INTERFACE;
  if(theString == MString("PARTICLE")) return (MInt)PARTICLE;
  if(theString == MString("SPECIES")) return (MInt)SPECIES;
  if(theString == MString("PATCH")) return (MInt)PATCH;
  if(theString == MString("CUTOFF")) return (MInt)CUTOFF;
  if(theString == MString("TOTALPRESSURE")) return (MInt)TOTALPRESSURE;
  if(theString == MString("DIVERGENCE")) return (MInt)DIVERGENCE;
  if(theString == MString("MEANSTRESS")) return (MInt)MEANSTRESS;
  if(theString == MString("SMOOTH")) return (MInt)SMOOTH;
  if(theString == MString("BAND")) return (MInt)BAND;

  // pascalm< RANS
  if(theString == MString("NORANS")) return (MInt)NORANS;
  if(theString == MString("RANS_SA")) return (MInt)RANS_SA;
  if(theString == MString("RANS_SA_DV")) return (MInt)RANS_SA_DV;
  if(theString == MString("RANS_FS")) return (MInt)RANS_FS;
  if(theString == MString("RANS_KOMEGA")) return (MInt)RANS_KOMEGA;
  if(theString == MString("RANS_SST")) return (MInt)RANS_SST;
  if(theString == MString("RANS_KEPSILON")) return (int)RANS_KEPSILON;
  //>pascalm RANS
  if(theString == MString("NAVIER_STOKES")) return (MInt)NAVIER_STOKES;
  if(theString == MString("EULER")) return (MInt)EULER;
  if(theString == MString("TINA_TC")) return (MInt)TINA_TC;
  if(theString == MString("SUZI_TC")) return (MInt)SUZI_TC;
  if(theString == MString("LOCD")) return (MInt)LOCD;
  if(theString == MString("HOCD")) return (MInt)HOCD;
  if(theString == MString("HOCD_LIMITED")) return (MInt)HOCD_LIMITED;
  if(theString == MString("HOCD_LIMITED_SLOPES")) return (MInt)HOCD_LIMITED_SLOPES;
  if(theString == MString("HOCD_LIMITED_SLOPES_MAN")) return (MInt)HOCD_LIMITED_SLOPES_MAN;
  if(theString == MString("AUSM")) return (MInt)AUSM;
  if(theString == MString("AUSMPLUS")) return (MInt)AUSMPLUS;
  if(theString == MString("SLAU")) return (MInt)SLAU;
  if(theString == MString("THREE_POINT")) return (MInt)THREE_POINT;
  if(theString == MString("FIVE_POINT")) return (MInt)FIVE_POINT;
  if(theString == MString("FIVE_POINT_STABILIZED")) return (MInt)FIVE_POINT_STABILIZED;
  if(theString == MString("FIVE_POINT_MULTISPECIES")) return (MInt)FIVE_POINT_MULTISPECIES;
  if(theString == MString("STANDARD")) return (MInt)STANDARD;
  if(theString == MString("LOCAL")) return (MInt)LOCAL;
  if(theString == MString("FAST")) return (MInt)FAST;
  if(theString == MString("CR1")) return (MInt)CR1;
  if(theString == MString("CR2")) return (MInt)CR2;
  if(theString == MString("HCR1")) return (MInt)HCR1;
  if(theString == MString("HCR2")) return (MInt)HCR2;
  if(theString == MString("HCR2_LIMITED")) return (MInt)HCR2_LIMITED;
  if(theString == MString("HCR2_FULLREINIT")) return (MInt)HCR2_FULLREINIT;
  if(theString == MString("SUS5CR")) return (MInt)SUS5CR;
  if(theString == MString("CR2PLUS")) return (MInt)CR2PLUS;
  if(theString == MString("RSU")) return (MInt)RSU;
  if(theString == MString("DL1")) return (MInt)DL1;
  if(theString == MString("DL2")) return (MInt)DL2;
  if(theString == MString("SUS_1")) return (MInt)SUS_1;
  if(theString == MString("SUS_1PLUS")) return (MInt)SUS_1PLUS;
  if(theString == MString("SUS_2")) return (MInt)SUS_2;
  if(theString == MString("SUS_WENO5")) return (MInt)SUS_WENO5;
  if(theString == MString("SUS_WENOO5PLUS")) return (MInt)SUS_WENO5PLUS;
  if(theString == MString("ELL")) return (MInt)ELL;
  if(theString == MString("no")) return (MInt)no;
  if(theString == MString("US1")) return (MInt)US1;
  if(theString == MString("UC3")) return (MInt)UC3;
  if(theString == MString("UC3_SB")) return (MInt)UC3_SB;
  if(theString == MString("UC5")) return (MInt)UC5;
  if(theString == MString("UC5_SB")) return (MInt)UC5_SB;
  if(theString == MString("UC11")) return (MInt)UC11;
  if(theString == MString("WENO5")) return (MInt)WENO5;
  if(theString == MString("WENO5_SB")) return (MInt)WENO5_SB;
  if(theString == MString("BACKWARDS_PAR")) return (MInt)BACKWARDS_PAR;
  if(theString == MString("ROTATING_LS")) return (MInt)ROTATING_LS;
  if(theString == MString("ROTATING_BNDRY")) return (MInt)ROTATING_BNDRY;
  if(theString == MString("BACKWARDS_INCREMENT")) return (MInt)BACKWARDS_INCREMENT;
  if(theString == MString("SYMMETRIC")) return (MInt)SYMMETRIC;
  if(theString == MString("PERIODIC")) return (MInt)PERIODIC;
  if(theString == MString("METHANE_2_STEP")) return (MInt)METHANE_2_STEP;
  if(theString == MString("METHANE_1_STEP")) return (MInt)METHANE_1_STEP;
  if(theString == MString("NONE")) return (MInt)NONE;
  if(theString == MString("MAIA_MAC_CORMACK")) return (MInt)MAIA_MAC_CORMACK;


  //#####################################################################
  //##      Postprocessing methods
  //#####################################################################

  if(theString == MString("PP_SLICE_AVERAGE")) return (MInt)PP_SLICE_AVERAGE;
  if(theString == MString("PP_REDUCE_TO_LEVEL_PRE")) return (MInt)PP_REDUCE_TO_LEVEL_PRE;
  if(theString == MString("PP_REDUCE_TO_LEVEL_POST")) return (MInt)PP_REDUCE_TO_LEVEL_POST;
  if(theString == MString("PP_REDUCE_TO_LEVEL_AVERAGES_PRE")) return (MInt)PP_REDUCE_TO_LEVEL_AVERAGES_PRE;
  if(theString == MString("PP_AVERAGE_PRE")) return (MInt)PP_AVERAGE_PRE;
  if(theString == MString("PP_AVERAGE_POST")) return (MInt)PP_AVERAGE_POST;
  if(theString == MString("PP_AVERAGE_IN")) return (MInt)PP_AVERAGE_IN;
  if(theString == MString("PP_COMPUTE_DIVERGENCEVELOCITY_PRE")) return (MInt)PP_COMPUTE_DIVERGENCEVELOCITY_PRE;
  if(theString == MString("PP_COMPUTE_DIVERGENCEVELOCITY_IN")) return (MInt)PP_COMPUTE_DIVERGENCEVELOCITY_IN;
  if(theString == MString("PP_COMPUTE_DIVERGENCEVELOCITY_POST")) return (MInt)PP_COMPUTE_DIVERGENCEVELOCITY_POST;
  if(theString == MString("PP_SPATIAL_AVERAGE_PRE")) return (MInt)PP_SPATIAL_AVERAGE_PRE;
  if(theString == MString("PP_SPATIAL_AVERAGE_POST")) return (MInt)PP_SPATIAL_AVERAGE_POST;
  if(theString == MString("PP_SPATIAL_AVERAGE_IN")) return (MInt)PP_SPATIAL_AVERAGE_IN;
  if(theString == MString("PP_MOVING_AVERAGE_PRE")) return (MInt)PP_MOVING_AVERAGE_PRE;
  if(theString == MString("PP_MOVING_AVERAGE_POST")) return (MInt)PP_MOVING_AVERAGE_POST;
  if(theString == MString("PP_MOVING_AVERAGE_IN")) return (MInt)PP_MOVING_AVERAGE_IN;
  if(theString == MString("PP_PROBE_POINT_PRE")) return (MInt)PP_PROBE_POINT_PRE;
  if(theString == MString("PP_PROBE_POINT_POST")) return (MInt)PP_PROBE_POINT_POST;
  if(theString == MString("PP_PROBE_POINT_IN")) return (MInt)PP_PROBE_POINT_IN;
  if(theString == MString("PP_PROBE_LINE_PRE")) return (MInt)PP_PROBE_LINE_PRE;
  if(theString == MString("PP_PROBE_LINE_POST")) return (MInt)PP_PROBE_LINE_POST;
  if(theString == MString("PP_PROBE_LINE_IN")) return (MInt)PP_PROBE_LINE_IN;
  if(theString == MString("PP_PROBE_LINE_PERIODIC_POST")) return (MInt)PP_PROBE_LINE_PERIODIC_POST;
  if(theString == MString("PP_PROBE_LINE_PERIODIC_IN")) return (MInt)PP_PROBE_LINE_PERIODIC_IN;
  if(theString == MString("PP_PROBE_ARB_LINE_PRE")) return (MInt)PP_PROBE_ARB_LINE_PRE;
  if(theString == MString("PP_PROBE_ARB_LINE_POST")) return (MInt)PP_PROBE_ARB_LINE_POST;
  if(theString == MString("PP_PROBE_ARB_LINE_IN")) return (MInt)PP_PROBE_ARB_LINE_IN;
  if(theString == MString("PP_PROBE_SLICE_PRE")) return (MInt)PP_PROBE_SLICE_PRE;
  if(theString == MString("PP_PROBE_SLICE_POST")) return (MInt)PP_PROBE_SLICE_POST;
  if(theString == MString("PP_PROBE_SLICE_IN")) return (MInt)PP_PROBE_SLICE_IN;
  if(theString == MString("PP_PROBE_ARB_SLICE_PRE")) return (MInt)PP_PROBE_ARB_SLICE_PRE;
  if(theString == MString("PP_PROBE_ARB_SLICE_POST")) return (MInt)PP_PROBE_ARB_SLICE_POST;
  if(theString == MString("PP_PROBE_ARB_SLICE_IN")) return (MInt)PP_PROBE_ARB_SLICE_IN;
  if(theString == MString("PP_WRITEPOINTS_IN")) return (MInt)PP_WRITEPOINTS_IN;
  if(theString == MString("PP_AVERAGE_SLICE_PRE")) return (MInt)PP_AVERAGE_SLICE_PRE;
  if(theString == MString("PP_TAUW_PRE")) return (MInt)PP_TAUW_PRE;
  if(theString == MString("PP_SUBTRACT_PERIODIC_FLUCTUATIONS")) return (MInt)PP_SUBTRACT_PERIODIC_FLUCTUATIONS;
  if(theString == MString("PP_SUBTRACT_MEAN")) return (MInt)PP_SUBTRACT_MEAN;
  if(theString == MString("PP_LOAD_AVERAGED_SOLUTION_PRE")) return (MInt)PP_LOAD_AVERAGED_SOLUTION_PRE;
  if(theString == MString("PP_COMPUTE_PRODUCTION_TERMS_PRE")) return (MInt)PP_COMPUTE_PRODUCTION_TERMS_PRE;
  if(theString == MString("PP_COMPUTE_DISSIPATION_TERMS_PRE")) return (MInt)PP_COMPUTE_DISSIPATION_TERMS_PRE;
  if(theString == MString("PP_WRITE_GRADIENTS")) return (MInt)PP_WRITE_GRADIENTS;
  if(theString == MString("PP_DECOMPOSE_CF")) return (MInt)PP_DECOMPOSE_CF;
  if(theString == MString("PP_SPRAY_STATS")) return (MInt)PP_SPRAY_STATS;
  if(theString == MString("PP_PARTICLE_SOLUTION")) return (MInt)PP_PARTICLE_SOLUTION;
  if(theString == MString("PP_POINT_SAMPLING_IN")) return (MInt)PP_POINT_SAMPLING_IN;
  if(theString == MString("PP_SURFACE_SAMPLING_IN")) return (MInt)PP_SURFACE_SAMPLING_IN;
  if(theString == MString("PP_VOLUME_SAMPLING_IN")) return (MInt)PP_VOLUME_SAMPLING_IN;
  if(theString == MString("PP_PARTICLE_STATISTICS")) return (MInt)PP_PARTICLE_STATISTICS;
  if(theString == MString("PP_ISO_TURBULENCE_STATISTICS")) return (MInt)PP_ISO_TURBULENCE_STATISTICS;
  if(theString == MString("PP_PL_ISO_TURBULENCE_STATISTICS")) return (MInt)PP_PL_ISO_TURBULENCE_STATISTICS;

  if(theString == MString("MINT")) return (MInt)MINT;
  if(theString == MString("MLONG")) return (MInt)MLONG;
  if(theString == MString("MFLOAT")) return (MInt)MFLOAT;
  //   if (theString == MString("DOUBLE"))
  //     return (MInt) DOUBLE;
  if(theString == MString("MSTRING")) return (MInt)MSTRING;
  if(theString == MString("MBOOL")) return (MInt)MBOOL;

  //##########################################################################################
  //##      LB ENUMS
  //##
  //##      contains:
  //##        - 1. LB discretization methods
  //##        - 2. LB computation methods
  //##        - 3. LB refinement methods
  //##        - 4. LB init methods
  //##        - 5. LB refilling schemes MB
  //##        - 6. LB interpolated bounce back schemes MB
  //##
  //##########################################################################################

  //#####################################################################
  //##      1. LB discretization methods
  //#####################################################################

  if(theString == MString("D2Q9")) return (MInt)D2Q9;
  if(theString == MString("D3Q15")) return (MInt)D3Q15;
  if(theString == MString("D3Q19")) return (MInt)D3Q19;
  if(theString == MString("D3Q27")) return (MInt)D3Q27;


  //#####################################################################
  //##      2. LB computation methods
  //#####################################################################

  if(theString == MString("MAIA_LATTICE_BGK")) return (MInt)MAIA_LATTICE_BGK;
  if(theString == MString("MAIA_LATTICE_BGK_INIT")) return (MInt)MAIA_LATTICE_BGK_INIT;
  if(theString == MString("MAIA_LATTICE_BGKI_SMAGORINSKY")) return (MInt)MAIA_LATTICE_BGKI_SMAGORINSKY;
  if(theString == MString("MAIA_LATTICE_BGKI_SMAGORINSKY2")) return (MInt)MAIA_LATTICE_BGKI_SMAGORINSKY2;
  if(theString == MString("MAIA_LATTICE_BGKI_SMAGO_WALL")) return (MInt)MAIA_LATTICE_BGKI_SMAGO_WALL;
  if(theString == MString("MAIA_LATTICE_BGKI_DYNAMIC_SMAGO")) return (MInt)MAIA_LATTICE_BGKI_DYNAMIC_SMAGO;
  if(theString == MString("MAIA_LATTICE_RBGK_DYNAMIC_SMAGO")) return (MInt)MAIA_LATTICE_RBGK_DYNAMIC_SMAGO;
  if(theString == MString("MAIA_LATTICE_BGKI_EULER_2D")) return (MInt)MAIA_LATTICE_BGKI_EULER_2D;
  if(theString == MString("MAIA_LATTICE_BGKC")) return (MInt)MAIA_LATTICE_BGKC;
  if(theString == MString("MAIA_LATTICE_RBGK")) return (MInt)MAIA_LATTICE_RBGK;
  if(theString == MString("MAIA_LATTICE_RBGK_SMAGORINSKY")) return (MInt)MAIA_LATTICE_RBGK_SMAGORINSKY;
  if(theString == MString("MAIA_LATTICE_MRT")) return (MInt)MAIA_LATTICE_MRT;
  if(theString == MString("MAIA_LATTICE_MRT2")) return (MInt)MAIA_LATTICE_MRT2;
  if(theString == MString("MAIA_LATTICE_CLB")) return (MInt)MAIA_LATTICE_CLB;
  if(theString == MString("MAIA_LATTICE_CLB_SMAGORINSKY")) return (MInt)MAIA_LATTICE_CLB_SMAGORINSKY;
  if(theString == MString("MAIA_LATTICE_MRT_SMAGORINSKY")) return (MInt)MAIA_LATTICE_MRT_SMAGORINSKY;
  if(theString == MString("MAIA_LATTICE_BGK_THERMAL")) return (MInt)MAIA_LATTICE_BGK_THERMAL;
  if(theString == MString("MAIA_LATTICE_BGK_INNERENERGY")) return (MInt)MAIA_LATTICE_BGK_INNERENERGY;
  if(theString == MString("MAIA_LATTICE_BGK_TOTALENERGY")) return (MInt)MAIA_LATTICE_BGK_TOTALENERGY;
  if(theString == MString("MAIA_LATTICE_BGK_TRANSPORT")) return (MInt)MAIA_LATTICE_BGK_TRANSPORT;
  if(theString == MString("MAIA_LATTICE_BGK_THERMAL_TRANSPORT")) return (MInt)MAIA_LATTICE_BGK_THERMAL_TRANSPORT;
  if(theString == MString("MAIA_LATTICE_BGK_INNERENERGY_TRANSPORT")) {
    return (MInt)MAIA_LATTICE_BGK_INNERENERGY_TRANSPORT;
  }
  if(theString == MString("MAIA_LATTICE_BGK_TOTALENERGY_TRANSPORT")) {
    return (MInt)MAIA_LATTICE_BGK_TOTALENERGY_TRANSPORT;
  }
  if(theString == MString("MAIA_LATTICE_CUMULANT")) return (MInt)MAIA_LATTICE_CUMULANT;
  if(theString == MString("MAIA_LATTICE_BGK_GUO_FORCING")) return (MInt)MAIA_LATTICE_BGK_GUO_FORCING;


  //#####################################################################
  //##      3. LB refinement methods
  //#####################################################################

  if(theString == MString("ROHDE")) return (MInt)ROHDE;
  if(theString == MString("FILIPPOVA")) return (MInt)FILIPPOVA;
  if(theString == MString("LINEAR_INTERPOLATION")) return (MInt)LINEAR_INTERPOLATION;
  if(theString == MString("QUADRATIC_INTERPOLATION")) return (MInt)QUADRATIC_INTERPOLATION;
  if(theString == MString("CUBIC_INTERPOLATION")) return (MInt)CUBIC_INTERPOLATION;

  //#####################################################################
  //##      4. LB init methods
  //#####################################################################

  // LB laminar inits (zero and direction-dependent)
  if(theString == MString("LB_FROM_ZERO_INIT")) return (MInt)LB_FROM_ZERO_INIT;
  if(theString == MString("LB_LAMINAR_INIT_PX")) return (MInt)LB_LAMINAR_INIT_PX;
  if(theString == MString("LB_LAMINAR_INIT_MX")) return (MInt)LB_LAMINAR_INIT_MX;
  if(theString == MString("LB_LAMINAR_INIT_PY")) return (MInt)LB_LAMINAR_INIT_PY;
  if(theString == MString("LB_LAMINAR_INIT_MY")) return (MInt)LB_LAMINAR_INIT_MY;
  if(theString == MString("LB_LAMINAR_INIT_PZ")) return (MInt)LB_LAMINAR_INIT_PZ;
  if(theString == MString("LB_LAMINAR_INIT_MZ")) return (MInt)LB_LAMINAR_INIT_MZ;

  // LB laminar inits (special geometries)
  if(theString == MString("LB_LAMINAR_CHANNEL_INIT")) return (MInt)LB_LAMINAR_CHANNEL_INIT;
  if(theString == MString("LB_LAMINAR_CYLINDER_INIT")) return (MInt)LB_LAMINAR_CYLINDER_INIT;
  if(theString == MString("LB_LAMINAR_PIPE_INIT")) return (MInt)LB_LAMINAR_PIPE_INIT;
  if(theString == MString("LB_TGV_INIT")) return (MInt)LB_TGV_INIT;
  if(theString == MString("LB_GAUSS_PULSE_INIT")) return (MInt)LB_GAUSS_PULSE_INIT;
  if(theString == MString("LB_GAUSS_DIFFUSION_INIT")) return (MInt)LB_GAUSS_DIFFUSION_INIT;
  if(theString == MString("LB_GAUSS_ADVECTION_INIT")) return (MInt)LB_GAUSS_ADVECTION_INIT;
  if(theString == MString("LB_SPINNING_VORTICIES_INIT")) return (MInt)LB_SPINNING_VORTICIES_INIT;
  if(theString == MString("LB_STEADY_VORTEX_INIT")) return (MInt)LB_STEADY_VORTEX_INIT;
  if(theString == MString("LB_CONVECTING_VORTEX_INIT")) return (MInt)LB_CONVECTING_VORTEX_INIT;

  // LB turbulent inits
  if(theString == MString("LB_TURBULENT_CHANNEL_INIT")) return (MInt)LB_TURBULENT_CHANNEL_INIT;
  if(theString == MString("LB_TURBULENT_MIXING_INIT")) return (MInt)LB_TURBULENT_MIXING_INIT;
  if(theString == MString("LB_TURBULENT_MIXING_FILTER_INIT")) return (MInt)LB_TURBULENT_MIXING_FILTER_INIT;
  if(theString == MString("LB_TURBULENT_PIPE_INIT")) return (MInt)LB_TURBULENT_PIPE_INIT;
  if(theString == MString("LB_TURBULENCE_ISOTROPIC_INIT")) return (MInt)LB_TURBULENCE_ISOTROPIC_INIT;

  // enum for init mehtods in lb adaptation
  if(theString == MString("INIT_COPYPASTE")) return (MInt)INIT_COPYPASTE;
  if(theString == MString("INIT_DUPUIS_FILIPPOVA")) return (MInt)INIT_DUPUIS_FILIPPOVA;

  if(theString == MString("POWERLAW")) return (MInt)POWERLAW;
  if(theString == MString("CARREAU")) return (MInt)CARREAU;

  //#####################################################################
  //##      5. LB Refilling Schemes MB
  //#####################################################################

  if(theString == MString("NORMAL_EXTRAPOLATION")) return (MInt)NORMAL_EXTRAPOLATION;
  if(theString == MString("AVERAGED_EXTRAPOLATION")) return (MInt)AVERAGED_EXTRAPOLATION;
  if(theString == MString("EQ_NEQ")) return (MInt)EQ_NEQ;
  if(theString == MString("VEL_CONSTRAINED_NORMAL_EXTRAPOLATION")) return (MInt)VEL_CONSTRAINED_NORMAL_EXTRAPOLATION;

  //#####################################################################
  //##      6. LB Interpolated Bounce Back Schemes MB
  //#####################################################################

  if(theString == MString("BOUZIDI_LINEAR")) return (MInt)BOUZIDI_LINEAR;
  if(theString == MString("BOUZIDI_QUADRATIC")) return (MInt)BOUZIDI_QUADRATIC;
  if(theString == MString("YU_LINEAR")) return (MInt)YU_LINEAR;
  if(theString == MString("YU_QUADRATIC")) return (MInt)YU_QUADRATIC;

  //##########################################################################################
  //##      DG ENUMS
  //##########################################################################################

  // Discontinuous Galerkin solver
  if(theString == MString("DG_POLY_LEGENDRE")) return (MInt)DG_POLY_LEGENDRE;
  if(theString == MString("DG_INTEGRATE_GAUSS")) return (MInt)DG_INTEGRATE_GAUSS;
  if(theString == MString("DG_INTEGRATE_GAUSS_LOBATTO")) return (MInt)DG_INTEGRATE_GAUSS_LOBATTO;
  if(theString == MString("DG_SYSEQN_EULER")) return (MInt)DG_SYSEQN_EULER;
  if(theString == MString("DG_SYSEQN_LINEARSCALARADV")) return (MInt)DG_SYSEQN_LINEARSCALARADV;
  if(theString == MString("DG_SYSEQN_ACOUSTICPERTURB")) return (MInt)DG_SYSEQN_ACOUSTICPERTURB;
  if(theString == MString("DG_ADAPTIVE_NONE")) return (MInt)DG_ADAPTIVE_NONE;
  if(theString == MString("DG_ADAPTIVE_TEST")) return (MInt)DG_ADAPTIVE_TEST;
  if(theString == MString("DG_ADAPTIVE_GRADIENT")) return (MInt)DG_ADAPTIVE_GRADIENT;
  // Time Integration Schemes
  if(theString == MString("DG_TIMEINTEGRATION_CARPENTER_4_5")) return (MInt)DG_TIMEINTEGRATION_CARPENTER_4_5;
  if(theString == MString("DG_TIMEINTEGRATION_TOULORGEC_4_8")) return (MInt)DG_TIMEINTEGRATION_TOULORGEC_4_8;
  if(theString == MString("DG_TIMEINTEGRATION_NIEGEMANN_4_14")) return (MInt)DG_TIMEINTEGRATION_NIEGEMANN_4_14;
  if(theString == MString("DG_TIMEINTEGRATION_NIEGEMANN_4_13")) return (MInt)DG_TIMEINTEGRATION_NIEGEMANN_4_13;
  if(theString == MString("DG_TIMEINTEGRATION_TOULORGEC_3_7")) return (MInt)DG_TIMEINTEGRATION_TOULORGEC_3_7;
  if(theString == MString("DG_TIMEINTEGRATION_TOULORGEF_4_8")) return (MInt)DG_TIMEINTEGRATION_TOULORGEF_4_8;

  if(theString == MString("BC_UNSET")) return (MInt)BC_UNSET;
  if(theString == MString("BC_DIRICHLET")) return (MInt)BC_DIRICHLET;
  if(theString == MString("BC_NEUMANN")) return (MInt)BC_NEUMANN;
  if(theString == MString("BC_ROBIN")) return (MInt)BC_ROBIN;
  if(theString == MString("BC_ISOTHERMAL")) return (MInt)BC_ISOTHERMAL;

  //##########################################################################################
  //##      FINITE CELL SOLVER ENUMS
  //##########################################################################################

  if(theString == MString("MAIA_LINEAR_ELASTIC")) return (MInt)MAIA_LINEAR_ELASTIC;
  if(theString == MString("MAIA_NONLINEAR_BAR")) return (MInt)MAIA_NONLINEAR_BAR;

  //##########################################################################################
  //##      STRUCTURED SOLVER ENUMS
  //##########################################################################################

  if(theString == MString("VENKATAKRISHNAN_MOD")) return (MInt)VENKATAKRISHNAN_MOD;
  if(theString == MString("VENKATAKRISHNAN")) return (MInt)VENKATAKRISHNAN;
  if(theString == MString("BARTH_JESPERSON")) return (MInt)BARTH_JESPERSON;
  if(theString == MString("MINMOD")) return (MInt)MINMOD;
  if(theString == MString("ALBADA")) return (MInt)ALBADA;
  if(theString == MString("PTHRC")) return (MInt)PTHRC;
  if(theString == MString("LEASTSQUARES")) return (MInt)LEASTSQUARES;

  // FV-Particle
  if(theString == MString("PART_EMITT_DIST_NONE")) return (MInt)PART_EMITT_DIST_NONE;
  if(theString == MString("PART_EMITT_DIST_UNIFORM")) return (MInt)PART_EMITT_DIST_UNIFORM;
  if(theString == MString("PART_EMITT_DIST_GAUSSIAN")) return (MInt)PART_EMITT_DIST_GAUSSIAN;

  if(theString == MString("SPRAY_ANGLE_MODEL_CONST")) return (MInt)SPRAY_ANGLE_MODEL_CONST;
  if(theString == MString("SPRAY_ANGLE_MODEL_HIROYASU_ARAI80")) return (MInt)SPRAY_ANGLE_MODEL_HIROYASU_ARAI80;
  if(theString == MString("SPRAY_ANGLE_MODEL_HIROYASU_ARAI90")) return (MInt)SPRAY_ANGLE_MODEL_HIROYASU_ARAI90;
  if(theString == MString("SPRAY_ANGLE_MODEL_BRACO_REITZ")) return (MInt)SPRAY_ANGLE_MODEL_BRACO_REITZ;
  if(theString == MString("SPRAY_ANGLE_MODEL_REITZ")) return (MInt)SPRAY_ANGLE_MODEL_REITZ;

  // solver surface type
  if(theString == MString("STL")) return (MInt)STL;
  if(theString == MString("ANALYTIC_BOX")) return (MInt)ANALYTIC_BOX;
  if(theString == MString("ANALYTIC_SPHERE")) return (MInt)ANALYTIC_SPHERE;

  // Phase description
  if(theString == MString("GAS")) return (MInt)GAS;
  if(theString == MString("LIQUID")) return (MInt)LIQUID;
  if(theString == MString("SOLID")) return (MInt)SOLID;

  // Injector types
  if(theString == MString("FULLCONE")) return (MInt)FULLCONE;
  if(theString == MString("HOLLOWCONE")) return (MInt)HOLLOWCONE;
  if(theString == MString("MULTIHOLE")) return (int)MULTIHOLE;
  if(theString == MString("MULTIHOLE_OPT")) return (int)MULTIHOLE_OPT;
  if(theString == MString("MULTIHOLE_TME")) return (int)MULTIHOLE_TME;


  // Dynamic load balancing
  if(theString == MString("DLB_PARTITION_DEFAULT")) return (MInt)DLB_PARTITION_DEFAULT;
  if(theString == MString("DLB_PARTITION_WEIGHT")) return (MInt)DLB_PARTITION_WEIGHT;
  if(theString == MString("DLB_PARTITION_SHIFT_OFFSETS")) return (MInt)DLB_PARTITION_SHIFT_OFFSETS;
  if(theString == MString("DLB_PARTITION_TEST")) return (MInt)DLB_PARTITION_TEST;

  // Acoustic extrapolation methods
  if(theString == MString("FWH")) return (MInt)FWH_METHOD;
  if(theString == MString("FWH_APE")) return (MInt)FWH_APE_METHOD;

  // Solver variable identifiers for sampling
  // FV Solver
  if(theString == MString("FV_SYSEQN_RANS")) return (MInt)FV_SYSEQN_RANS;
  if(theString == MString("FV_SYSEQN_EEGAS")) return (MInt)FV_SYSEQN_EEGAS;
  if(theString == MString("FV_SYSEQN_NS")) return (MInt)FV_SYSEQN_NS;
  if(theString == MString("FV_SYSEQN_DETCHEM")) return (MInt)FV_SYSEQN_DETCHEM;
  if(theString == MString("FV_PV")) return (MInt)FV_PV;
  if(theString == MString("FV_VORT")) return (MInt)FV_VORT;
  if(theString == MString("FV_HEAT_RELEASE")) return (MInt)FV_HEAT_RELEASE;
  // DG Solver
  if(theString == MString("DG_VARS")) return (MInt)DG_VARS;
  if(theString == MString("DG_NODEVARS")) return (MInt)DG_NODEVARS;
  if(theString == MString("DG_SOURCETERMS")) return (MInt)DG_SOURCETERMS;

  if(theString == MString("LB_PV")) return (MInt)LB_PV;
  // Unified run loop - coupler types
  if(theString == MString("COUPLER_LS_FV_MB")) return (MInt)COUPLER_LS_FV_MB;
  if(theString == MString("COUPLER_LS_FV")) return (MInt)COUPLER_LS_FV;
  if(theString == MString("COUPLER_FV_MULTILEVEL")) return (MInt)COUPLER_FV_MULTILEVEL;
  if(theString == MString("COUPLER_FV_ZONAL_RTV")) return (MInt)COUPLER_FV_ZONAL_RTV;
  if(theString == MString("COUPLER_FV_MULTILEVEL_INTERPOLATION")) return (MInt)COUPLER_FV_MULTILEVEL_INTERPOLATION;
  if(theString == MString("COUPLER_FV_ZONAL_STG")) return (MInt)COUPLER_FV_ZONAL_STG;
  if(theString == MString("COUPLER_CARTESIAN_INTERPOLATION")) return (MInt)COUPLER_CARTESIAN_INTERPOLATION;
  if(theString == MString("COUPLER_FV_DG_APE")) return (MInt)COUPLER_FV_DG_APE;
  if(theString == MString("COUPLER_LS_FV_COMBUSTION")) return (MInt)COUPLER_LS_FV_COMBUSTION;
  if(theString == MString("COUPLER_LS_LB")) return (MInt)COUPLER_LS_LB;
  if(theString == MString("COUPLER_LS_LB_PARTICLE")) return (MInt)COUPLER_LS_LB_PARTICLE;
  if(theString == MString("COUPLER_LB_LPT")) return (MInt)COUPLER_LB_LPT;
  if(theString == MString("COUPLER_FV_PARTICLE")) return (MInt)COUPLER_FV_PARTICLE;
  if(theString == MString("COUPLER_LS_LB")) return (MInt)COUPLER_LS_LB;
  if(theString == MString("COUPLER_LS_LB_SURFACE")) return (MInt)COUPLER_LS_LB_SURFACE;
  if(theString == MString("COUPLER_LS_FV_COMBUSTION")) return (MInt)COUPLER_LS_FV_COMBUSTION;
  if(theString == MString("COUPLER_FV_MB_ZONAL")) return (MInt)COUPLER_FV_MB_ZONAL;
  if(theString == MString("COUPLER_LB_FV_EE_MULTIPHASE")) return (MInt)COUPLER_LB_FV_EE_MULTIPHASE;
  if(theString == MString("COUPLER_LB_LB")) return (MInt)COUPLER_LB_LB;
  if(theString == MString("COUPLER_LB_DG_APE")) return (MInt)COUPLER_LB_DG_APE;
  if(theString == MString("COUPLER_LB_RB")) return (MInt)COUPLER_LB_RB;

  // Unified run loop - postprocessing types
  if(theString == MString("POSTPROCESSING_FV")) return (MInt)POSTPROCESSING_FV;
  if(theString == MString("POSTPROCESSING_LS")) return (MInt)POSTPROCESSING_LS;
  if(theString == MString("POSTPROCESSING_DG")) return (MInt)POSTPROCESSING_DG;
  if(theString == MString("POSTPROCESSING_LB")) return (MInt)POSTPROCESSING_LB;
  if(theString == MString("POSTPROCESSING_FVLPT")) return (MInt)POSTPROCESSING_FVLPT;
  if(theString == MString("POSTPROCESSING_LBLPT")) return (MInt)POSTPROCESSING_LBLPT;

  // Unified run loop - recipe types
  if(theString == MString("RECIPE_INTRASTEP")) return (MInt)RECIPE_INTRASTEP;
  if(theString == MString("RECIPE_BASE")) return (MInt)RECIPE_BASE;
  if(theString == MString("RECIPE_ITERATION")) return (MInt)RECIPE_ITERATION;


  // Detailed chemistry transport models
  if(theString == MString("Multi")) return (MInt)Multi;
  if(theString == MString("Mix")) return (MInt)Mix;

  // Viscosity laws
  if(theString == MString("SUTHERLAND")) return (MInt)SUTHERLAND;
  if(theString == MString("CONSTANT")) return (MInt)CONSTANT;

  if(theString == MString("PARTICLE_COUNT")) return (MInt)PARTICLE_COUNT;
  if(theString == MString("PARTICLE_FLOAT")) return (MInt)PARTICLE_FLOAT;
  if(theString == MString("PARTICLE_INT")) return (MInt)PARTICLE_INT;
  if(theString == MString("SOURCE_TERMS")) return (MInt)SOURCE_TERMS;
  if(theString == MString("FLOW_FIELD")) {
    return (MInt)FLOW_FIELD;
  }
  if(theString == MString("CHECK_ADAP")) {
    return (MInt)CHECK_ADAP;
  }
  if(theString == MString("VELOCITY_SLOPES")) {
    return (MInt)VELOCITY_SLOPES;
  }

  // Finite cell interpolation methods
  if(theString == MString("EQUIDIST_LAGRANGE_INTERP")) return (MInt)EQUIDIST_LAGRANGE_INTERP;
  if(theString == MString("LAGRANGE_INTERP")) return (MInt)LAGRANGE_INTERP;
  if(theString == MString("EQUIDIST_LEGENDRE_INTERP")) return (MInt)EQUIDIST_LEGENDRE_INTERP;
  if(theString == MString("LEGENDRE_INTERP")) return (MInt)LEGENDRE_INTERP;

  // FW-H methods in ACA solver
  if(theString == MString("MAIA_FWH_FREQUENCY")) return (MInt)MAIA_FWH_FREQUENCY;
  if(theString == MString("MAIA_FWH_TIME")) return (MInt)MAIA_FWH_TIME;

  cerr << "In function string2enum(): No enum definition for string '" << theString << "' found !" << endl;
  return -1;
}
