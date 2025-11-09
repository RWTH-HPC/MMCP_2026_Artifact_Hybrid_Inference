// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <doctest/doctest.h>

#include "LB/lbfunctions.h"

// Combine lattice properties and physical model in one type
template <MInt nDim_, MInt nDist_, MBool compressible_>
struct LbModel {
  static constexpr MInt nDim = nDim_;
  static constexpr MInt nDist = nDist_;
  static constexpr MBool compressible = compressible_;
};

TEST_SUITE_BEGIN("LB_FUNCTIONS");

// Check if macro vars -> eq dist -> macro vars is an identity
TEST_CASE_TEMPLATE_DEFINE("LB equlibrium function identity", Model, test_id_descriptor) {
  constexpr MInt nDim = Model::nDim;
  constexpr MInt nDist = Model::nDist;
  constexpr MBool compressible = Model::compressible;

  // Input
  MFloat zerothMoment{};
  std::array<MFloat, nDim> firstMoment{};

  // Define different input data as sub cases
  // Calculations/Checks will be performed for each
  SUBCASE("First moment zero") {
    zerothMoment = 1.0;
    firstMoment.fill(0.0);
  }

  SUBCASE("First moment unity") {
    zerothMoment = 1.0;
    firstMoment.fill(1.0);
  }

  // Transform to equilibrium distribution and back..
  std::array<MFloat, nDist> dist{};
  lbfunc::calcEqDists<nDim, nDist, compressible>(zerothMoment, firstMoment.data(), dist.data());
  MFloat zerothMoment_res{};
  std::array<MFloat, nDim> firstMoment_res{};
  lbfunc::calcMacroVars<nDim, nDist, compressible>(dist.data(), zerothMoment_res, firstMoment_res.data());

  CHECK(zerothMoment_res == doctest::Approx(zerothMoment));

  for(MInt i = 0; i < nDim; i++) {
    CHECK(firstMoment_res[i] == doctest::Approx(firstMoment[i]));
  }
}

// Compressible cases
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbModel<2, 9, true>);
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbModel<3, 19, true>);
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbModel<3, 27, true>);

// Incompressible cases
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbModel<2, 9, false>);
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbModel<3, 19, false>);
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbModel<3, 27, false>);

TEST_SUITE_END();
