// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Each unit test has to include the following
#include <doctest/doctest.h>

// Now, simple include everything needed for the following tests
#include "LB/lblatticedescriptor.h"

// Test suite may be set-up (not necessary). It is useful if you would like to
// call an isolated set of test cases by calling
//    ./test_maia --test-suite=LB_LATTICEDESCRIPTOR
TEST_SUITE_BEGIN("LB_LATTICEDESCRIPTOR");

// Usual test case with no sub cases
// First, create the test case with a certain description
TEST_CASE("Test dirFld2 for duplicated values") {
  // Create/initialize the stuff you need for your test case/s
  constexpr MInt dim = 3;
  MBool check = true;
  for(MInt i = 0; i < dim * dim; i++) {
    const MInt ii = i / dim;
    const MInt ij = i % dim;
    const MInt a = lbDescriptor::dirFld2[ii][ij];
    for(MInt j = i + 1; j < dim * dim; j++) {
      const MInt ji = j / dim;
      const MInt jj = j % dim;
      const MInt b = lbDescriptor::dirFld2[ji][jj];
      check = check && (a != b);
    }
  }
  // Check the relevant stuff using e.g. CHECK
  CHECK(check);

  /*
   * Further operation/calculation could be done here, too.
   * These operations would then only affect the CHECKs following afterwards.
   */

  /*
   * Additional sub cases, i.e., 'branches' of test can be introduced using
   * SUBCASE("descripton subcase, i.e., branch of testing"){
   *    // Further steps/calculations
   *    // Further CHECKs
   * }
   */
}

// Test case template definition using two SUBCASES
TEST_CASE_TEMPLATE_DEFINE("LB lattice descriptor", Ld, test_id_descriptor) {
  // Everything here, would affect (be usable in) both SUBCASEs !
  SUBCASE("lastId") { CHECK(Ld::lastId() == Ld::q() - 1); }
  SUBCASE("oppositeDist") {
    for(MInt i = 0; i < Ld::q() / 2; i++) {
      const MInt iOpp = Ld::oppositeDist(i);
      const MInt iOppOpp = Ld::oppositeDist(iOpp);
      CHECK(i != iOpp);
      CHECK(i == iOppOpp);
    }
  }
}
// Now, we use the template test case, whose id is stored in test_id_descriptor
// to instantiate it multiple times with different class types
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbLatticeDescriptor<2, 9>);
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbLatticeDescriptor<3, 19>);
TEST_CASE_TEMPLATE_INVOKE(test_id_descriptor, LbLatticeDescriptor<3, 27>);

// We end this suite, further test and/or suite may follow
TEST_SUITE_END();
