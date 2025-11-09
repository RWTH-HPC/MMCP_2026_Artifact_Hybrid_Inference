// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_VARIABLES_H_
#define MAIA_VARIABLES_H_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file defines the index classes for the access of within
/// MAIA's arrays of variables (e.g. Cell.m_variables, solver.m_variables,...)
////////////////////////////////////////////////////////////////////////////////
#include <array>
#include <limits>
#include "INCLUDE/maiatypes.h"
////////////////////////////////////////////////////////////////////////////////

/*!\class MConservativeVariables
   \date begin: 00.05.19 change 00.05.19
   \\author changed by Pascal Meysonnat to deal also with RANS

   \brief Storage of the Position of the Conservative Variables (RHO, RHO_VV, RHO_E)
   in the value vectors of the solvers ans surfaces

   If you search for a value in a value vector you only have to type in the name
   of the value you are searching (e.g. variables[RHO_E])
   \warning Be careful! Some value vectors include the conservative, others the
   primitive Variables
*/
template <MInt nDim>
class MConservativeVariables {
 public:
  //! Sets the position for the conservative variables
  explicit MConservativeVariables(MInt noSpecies, MInt noRANSEq = 0) {
    m_noSpecies = noSpecies;
    m_noRansEquations = noRANSEq;
    for(MInt i = 0; i < nDim; i++) {
      RHO_VV[i] = RHO_U + i;
    }
    RHO_Y = new MInt[m_noSpecies];
    for(MInt i = 0; i < m_noSpecies; i++) {
      RHO_Y[i] = RHO_Z + i;
    }
    RANS_FIRST = nDim + 2;
    if(m_noSpecies != 0) {
      RANS_FIRST = nDim + 2 + m_noSpecies;
    }
    if(m_noRansEquations > 0) {
      RANS_VAR = new MInt[m_noRansEquations];
      for(MInt i = 0; i < m_noRansEquations; ++i) {
        RANS_VAR[i] = RANS_FIRST + i;
      }
    } else {
      RANS_VAR = nullptr;
    }
    noVariables = nDim + 2 + m_noSpecies + m_noRansEquations;
  }

  ~MConservativeVariables() {
    delete[] RHO_Y;
    if(RANS_VAR != nullptr) {
      delete[] RANS_VAR;
    }
  }

  //! Position of RHO_U
  static constexpr MInt RHO_U = 0;
  //! Position of RHO_V
  static constexpr MInt RHO_V = 1;
  //! Position of RHO_W
  static constexpr MInt RHO_W = 2;
  //! Pointer for the velocities so you can use them in a loop
  std::array<MInt, nDim> RHO_VV;
  //! Position of RHO_E
  static constexpr MInt RHO_E = nDim;
  //! Position of RHO
  static constexpr MInt RHO = nDim + 1;
  //! Position of RHO_Z
  static constexpr MInt RHO_Z = nDim + 2;
  //! Position of RHO_C
  static constexpr MInt RHO_C = nDim + 2;
  //! Position of RHO_Yi
  MInt* RHO_Y = nullptr;
  //! first Position of RANS Variables
  MInt RANS_FIRST;
  //! Position of RANS Variables
  MInt* RANS_VAR = nullptr;


  //! The Nr. of Conservative Variables (nDim + 2)
  MInt noVariables;
  MInt m_noSpecies;
  MInt m_noRansEquations;
  MFloat rhoUInfinity{}, rhoVInfinity{}, rhoWInfinity{}, rhoEInfinity{}, rhoInfinity{};
  std::array<MFloat, nDim> ransInfinity{};
  std::array<MFloat, nDim> rhoVVInfinity{};
};

/*!\class MPrimitiveVariables
   \date begin: 00.05.19 change 00.05.19
   \brief Storage of the Position of the Primitive Variables (u, v, w, T, p) in
   the value vectors of the solvers ans surfaces

   If you search for a value in a value vector you only have to type in the name
   of the value you are searching (e.g. variables[T])
   \warning Be careful! Some value vectors include the conservative, others the
   primitive Variables
*/
template <MInt nDim>
class MPrimitiveVariables {
 public:
  //! Sets the position for the primitive variables
  explicit MPrimitiveVariables(MInt noSpecies, MInt noRANSEq = 0) {
    m_noSpecies = noSpecies;
    m_noRansEquations = noRANSEq;

    Y = new MInt[m_noSpecies];
    for(MInt i = 0; i < m_noSpecies; i++) {
      Y[i] = Z + i;
    }
    // G = nDim+1;

    RANS_FIRST = nDim + 2;
    if(m_noSpecies != 0) {
      RANS_FIRST = nDim + 2 + m_noSpecies;
    }
    if(m_noRansEquations > 0) {
      RANS_VAR = new MInt[m_noRansEquations];
      for(MInt i = 0; i < m_noRansEquations; ++i) {
        RANS_VAR[i] = RANS_FIRST + i;
      }
    } else {
      RANS_VAR = nullptr;
    }

    noVariables = nDim + 2 + m_noSpecies + m_noRansEquations;
  }
  ~MPrimitiveVariables() {
    delete[] Y;
    if(RANS_VAR != nullptr) {
      delete[] RANS_VAR;
    }
  }

  //! Position of U
  static constexpr MInt U = 0;
  //! Position of V
  static constexpr MInt V = 1;
  //! Position of W
  static constexpr MInt W = 2;
  //! Pointer for the velocities so you can use them in a loop
#if !defined(MAIA_PGI_COMPILER) && defined(NDEBUG) && !defined(MAIA_SANITIZE_ADDRESS) && !defined(MAIA_INTEL_COMPILER)
  static constexpr std::array<MInt, 3> VV = {0, 1, 2};
#else
  std::array<MInt, 3> VV = {0, 1, 2};
#endif
  //! Position of P
  static constexpr MInt P = nDim + 1;
  //! Position of RHO (equal to P in this case)
  static constexpr MInt RHO = nDim;
  //! Position of T
  static constexpr MInt T = nDim + 1;
  //! Position of Z
  static constexpr MInt Z = nDim + 2;
  //! Position of C
  static constexpr MInt C = nDim + 2;
  //! Position of Yi
  MInt* Y;
  //! first Position of RANS Variables
  MInt RANS_FIRST;
  //! Position of RANS Variables
  MInt* RANS_VAR;

  //! The Nr. of primitive variables (nDim + 2)
  MInt noVariables;
  MInt m_noSpecies;
  MInt m_noRansEquations;
  MFloat UInfinity{}, VInfinity{}, WInfinity{}, PInfinity{}, TInfinity{};
  std::array<MFloat, nDim> VVInfinity;
  std::array<MFloat, nDim> ransInfinity;
  MFloat DthInfinity{}, muInfinity{}, DInfinity{};
};

////////////////////////////////////////////////////////////////////////////////
/// Classes with constant number of space dimensions nd

namespace maia {

namespace fv {

namespace variables {
/// Error value: produces a segmentation fault if used as index.
static const MInt Segfault = std::numeric_limits<MInt>::min();
} // namespace variables

/// \brief Static indices for accessing conservative variables
/// in nd spatial dimensions
template <MInt nd>
struct ConservativeVariables {
  static constexpr MInt RHO_U = 0;
  static constexpr MInt RHO_V = 1;
  static constexpr MInt RHO_W = nd == 3 ? 2 : variables::Segfault;
  static constexpr std::array<MInt, 3> RHO_VV = {0, 1, 2};
  static constexpr MInt RHO_E = nd;
  static constexpr MInt RHO = nd + 1;
  const MInt RHO_Z;
  const MInt RHO_C;
  static constexpr MInt RHO_N = nd + 2;

  const MInt m_noSpecies;
  const MInt m_noRansEquations;
  const MInt noVariables;

  MInt* RHO_NN = nullptr;
  MInt* RHO_Y = nullptr;

  explicit ConservativeVariables(const MInt noSpecies, const MInt noRans);
  ~ConservativeVariables();
};

/// \brief Static indices for accessing primitive variables
/// in nd spatial dimensions
template <MInt nd>
struct PrimitiveVariables {
  static constexpr MInt U = 0;
  static constexpr MInt V = 1;
  static constexpr MInt W = nd == 3 ? 2 : variables::Segfault;
  static constexpr std::array<MInt, 3> VV = {0, 1, 2};
  static constexpr MInt RHO = nd;
  static constexpr MInt P = nd + 1;
  static constexpr MInt T = nd + 1;

  const MInt Z;
  const MInt C;
  static constexpr MInt N = nd + 2;

  const MInt m_noSpecies;
  const MInt m_noRansEquations;
  const MInt noVariables;

  MInt* NN = nullptr;
  MInt* Y = nullptr;

  explicit PrimitiveVariables(const MInt noSpecies, const MInt noRans);
  ~PrimitiveVariables();
};

} // namespace fv
} // namespace maia
#endif
