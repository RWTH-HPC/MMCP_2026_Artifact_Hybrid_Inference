// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBLATTICEDESCRIPTOR_H
#define LBLATTICEDESCRIPTOR_H

#include "INCLUDE/maiaconstants.h"
#include "UTIL/functions.h"

// Description: This file contains arrays depending on number of space
// dimensions D and of populations Q. All these template arrays are wrapped in
// the struct LbLatticeDescriptorBase if only depending on D and in
// LbLatticeDescriptor if depending and D and Q.
// Usage: Defining a template specific typedef once in a class all arrays can be
// easily accesed, p.e.:
//    using Ld = LbLatticeDescriptor<nDim, nDist>;
//    Ld::idFld(i, j);

//--declarations of lattice descriptor------------------------------------------
namespace lbDescriptor {

//-- Dx --
/// Convert directions into distribution id
inline constexpr MInt dirFld2[3][3] = {{6, 0, 7}, {2, 8, 3}, {5, 1, 4}};
inline constexpr MInt dirFld3[3][3][3] = {{{18, 6, 19}, {10, 0, 11}, {20, 7, 21}},
                                          {{14, 2, 15}, {4, 26, 5}, {16, 3, 17}},
                                          {{22, 8, 23}, {12, 1, 13}, {24, 9, 25}}};

/// Basic movements for the distributions , 0=backwards, 1=stop, 2=forwards
template <MInt D>
inline constexpr MInt idFld[POWX(3, D)][D] = {};

/// Directions of the interpolation neighbors [positionChild][positionIntepolationParent]
template <MInt D>
inline constexpr MInt intNghbrArray[POWX(2, D)][POWX(2, D)] = {};

/// Linear interpolation coefficients
template <MInt D>
inline constexpr MFloat linearInterpolationCoefficients[POWX(2, D)][POWX(2, D)] = {};

/// Ids for the calculation of the first term of the maxwellian
template <MInt D>
inline constexpr MInt mFld1[24] = {};
/// Ids for the calculation of the second term of the maxwellian
inline constexpr MInt mFld2[24] = {0, 2, 4, 0, 2, 5, 0, 3, 4, 0, 3, 5, 1, 2, 4, 1, 2, 5, 1, 3, 4, 1, 3, 5};

/// Distributions with components in negative axis direction [axis][distribution]
template <MInt D>
inline constexpr MInt nFld[D][POWX(3, D - 1)] = {};
/// Distributions with components in positive axis direction [axis][distribution]
template <MInt D>
inline constexpr MInt pFld[D][POWX(3, D - 1)] = {};

/// Distribution ids required to form a vector for nodal connectivity
template <MInt D>
inline constexpr MInt nodalConnectivityVector[D * POWX(2, D - 1)][2] = {};

/// Opposite distribution
template <MInt D>
inline constexpr MInt oppositeDist[POWX(3, D)] = {};

/// Direction of motion of each PPDF, needed for bounce-back schemes
template <MInt D>
inline constexpr MFloat ppdfDir[POWX(3, D)][D] = {};

//-- DxQy --
/// Number of distributions - (1d,2d,3d)
template <MInt D, MInt Q>
inline constexpr MInt distFld[3] = {};

///
template <MInt D, MInt Q>
inline constexpr MInt distType[Q] = {};

/// Number of distributions which have one space direction in common
template <MInt D, MInt Q>
inline constexpr MInt dxQyFld = {};

/// taylor-polynom coefficients for equilibrium calculation (0:rest, 1:face, 2:edge, 3:corner)
template <MInt D, MInt Q>
inline constexpr MFloat tp[4] = {};
} // namespace lbDescriptor

//--definitions of lattice descriptor for specific dxqy-------------------------
namespace lbDescriptor {
//-- D2Qx --
template <>
inline constexpr MInt idFld<2>[9][2] = {{0, 1}, {2, 1}, {1, 0}, {1, 2}, {2, 2}, {2, 0}, {0, 0}, {0, 2}, {1, 1}};
template <>
inline constexpr MInt intNghbrArray<2>[4][4] = {{dirFld2[0][0], dirFld2[1][0], dirFld2[0][1], dirFld2[1][1]}, //
                                                {dirFld2[1][0], dirFld2[2][0], dirFld2[1][1], dirFld2[2][1]}, //
                                                {dirFld2[0][1], dirFld2[1][1], dirFld2[0][2], dirFld2[1][2]}, //
                                                {dirFld2[1][1], dirFld2[2][1], dirFld2[1][2], dirFld2[2][2]}};
template <>
inline constexpr MFloat linearInterpolationCoefficients<2>[4][4] = {{F1B16, F3B16, F3B16, F9B16},
                                                                    {F3B16, F1B16, F9B16, F3B16},
                                                                    {F3B16, F9B16, F1B16, F3B16},
                                                                    {F9B16, F3B16, F3B16, F1B16}};
template <>
inline constexpr MInt mFld1<2>[24] = {1,  3,  1,  2,  0,  2,  0,  3,  -1, -1, -1, -1,
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
template <>
inline constexpr MInt nFld<2>[2][3] = {{7, 6, 0}, {6, 2, 5}};
template <>
inline constexpr MInt pFld<2>[2][3] = {{1, 4, 5}, {7, 3, 4}};
template <>
inline constexpr MInt nodalConnectivityVector<2>[4][2] = {{6, 7},  // left edge
                                                          {4, 5},  // right edge
                                                          {5, 6},  // bottom edge
                                                          {7, 4}}; // top edge
template <>
inline constexpr MInt oppositeDist<2>[9] = {1, 0, 3, 2, 6, 7, 4, 5, 8};
template <>
inline constexpr MFloat ppdfDir<2>[9][2] = {{-1.0, 0.0}, {1.0, 0.0},   {0.0, -1.0}, {0.0, 1.0}, {1.0, 1.0},
                                            {1.0, -1.0}, {-1.0, -1.0}, {-1.0, 1.0}, {0.0, 0.0}};

//-- D3Qx --
template <>
inline constexpr MInt idFld<3>[27][3] = {{0, 1, 1}, {2, 1, 1}, // Dist. 0, 1
                                         {1, 0, 1}, {1, 2, 1}, // Dist. 2, 3
                                         {1, 1, 0}, {1, 1, 2}, // Dist. 4, 5


                                         {0, 0, 1}, {0, 2, 1}, // Dist. 6, 7
                                         {2, 0, 1}, {2, 2, 1}, // Dist. 8, 9

                                         {0, 1, 0}, {0, 1, 2}, // Dist. 10, 11
                                         {2, 1, 0}, {2, 1, 2}, // Dist. 12, 13

                                         {1, 0, 0}, {1, 0, 2}, // Dist. 14, 15
                                         {1, 2, 0}, {1, 2, 2}, // Dist. 16, 17

                                         {0, 0, 0}, {0, 0, 2}, // Dist. 18, 19
                                         {0, 2, 0}, {0, 2, 2}, // Dist. 20, 21
                                         {2, 0, 0}, {2, 0, 2}, // Dist. 22, 23
                                         {2, 2, 0}, {2, 2, 2}, // Dist. 24, 25

                                         {1, 1, 1}}; // Dist. 26
template <>
inline constexpr MInt intNghbrArray<3>[8][8] = {
    {dirFld3[0][0][0], dirFld3[1][0][0], dirFld3[0][1][0], dirFld3[1][1][0], dirFld3[0][0][1], dirFld3[1][0][1],
     dirFld3[0][1][1], dirFld3[1][1][1]},
    {dirFld3[1][0][0], dirFld3[2][0][0], dirFld3[1][1][0], dirFld3[2][1][0], dirFld3[1][0][1], dirFld3[2][0][1],
     dirFld3[1][1][1], dirFld3[2][1][1]},
    {dirFld3[0][1][0], dirFld3[1][1][0], dirFld3[0][2][0], dirFld3[1][2][0], dirFld3[0][1][1], dirFld3[1][1][1],
     dirFld3[0][2][1], dirFld3[1][2][1]},
    {dirFld3[1][1][0], dirFld3[2][1][0], dirFld3[1][2][0], dirFld3[2][2][0], dirFld3[1][1][1], dirFld3[2][1][1],
     dirFld3[1][2][1], dirFld3[2][2][1]},
    {dirFld3[0][0][1], dirFld3[1][0][1], dirFld3[0][1][1], dirFld3[1][1][1], dirFld3[0][0][2], dirFld3[1][0][2],
     dirFld3[0][1][2], dirFld3[1][1][2]},
    {dirFld3[1][0][1], dirFld3[2][0][1], dirFld3[1][1][1], dirFld3[2][1][1], dirFld3[1][0][2], dirFld3[2][0][2],
     dirFld3[1][1][2], dirFld3[2][1][2]},
    {dirFld3[0][1][1], dirFld3[1][1][1], dirFld3[0][2][1], dirFld3[1][2][1], dirFld3[0][1][2], dirFld3[1][1][2],
     dirFld3[0][2][2], dirFld3[1][2][2]},
    {dirFld3[1][1][1], dirFld3[2][1][1], dirFld3[1][2][1], dirFld3[2][2][1], dirFld3[1][1][2], dirFld3[2][1][2],
     dirFld3[1][2][2], dirFld3[2][2][2]}};
template <>
inline constexpr MFloat linearInterpolationCoefficients<3>[8][8] = {
    {F1B64, F3B64, F3B64, F9B64, F3B64, F9B64, F9B64, F27B64},
    {F3B64, F1B64, F9B64, F3B64, F9B64, F3B64, F27B64, F9B64},
    {F3B64, F9B64, F1B64, F3B64, F9B64, F27B64, F3B64, F9B64},
    {F9B64, F3B64, F3B64, F1B64, F27B64, F9B64, F9B64, F3B64},
    {F3B64, F9B64, F9B64, F27B64, F1B64, F3B64, F3B64, F9B64},
    {F9B64, F3B64, F27B64, F9B64, F3B64, F1B64, F9B64, F3B64},
    {F9B64, F27B64, F3B64, F9B64, F3B64, F9B64, F1B64, F3B64},
    {F27B64, F9B64, F9B64, F3B64, F9B64, F3B64, F3B64, F1B64}};
template <>
inline constexpr MInt mFld1<3>[24] = {0, 2, 0, 3, 1, 2, 1, 3, 0, 4, 0, 5, 1, 4, 1, 5, 2, 4, 2, 5, 3, 4, 3, 5};
template <>
inline constexpr MInt nFld<3>[3][9] = {
    {0, 6, 7, 10, 11, 18, 19, 20, 21}, {2, 6, 8, 14, 15, 18, 19, 22, 23}, {4, 10, 12, 14, 16, 18, 20, 22, 24}};
template <>
inline constexpr MInt pFld<3>[3][9] = {
    {1, 8, 9, 12, 13, 22, 23, 24, 25}, {3, 7, 9, 16, 17, 20, 21, 24, 25}, {5, 11, 13, 15, 17, 19, 21, 23, 25}};
template <>
inline constexpr MInt nodalConnectivityVector<3>[12][2] = {{18, 19},  // edge of Dist. 6
                                                           {20, 21},  // edge of Dist. 7
                                                           {22, 23},  // edge of Dist. 8
                                                           {24, 25},  // edge of Dist. 9
                                                           {18, 20},  // edge of Dist. 10
                                                           {19, 21},  // edge of Dist. 11
                                                           {22, 24},  // edge of Dist. 12
                                                           {23, 25},  // edge of Dist. 13
                                                           {18, 22},  // edge of Dist. 14
                                                           {19, 23},  // edge of Dist. 15
                                                           {20, 24},  // edge of Dist. 16
                                                           {21, 25}}; // edge of Dist. 17

template <>
inline constexpr MInt oppositeDist<3>[27] = {1,  0,  3,  2,  5,  4,  9,  8,  7,  6,  13, 12, 11, 10,
                                             17, 16, 15, 14, 25, 24, 23, 22, 21, 20, 19, 18, 26};
template <>
inline constexpr MFloat ppdfDir<3>[27][3] = {
    {-1.0, 0.0, 0.0},   {1.0, 0.0, 0.0},   {0.0, -1.0, 0.0},  {0.0, 1.0, 0.0},  {0.0, 0.0, -1.0},  {0.0, 0.0, 1.0},
    {-1.0, -1.0, 0.0},  {-1.0, 1.0, 0.0},  {1.0, -1.0, 0.0},  {1.0, 1.0, 0.0},  {-1.0, 0.0, -1.0}, {-1.0, 0.0, 1.0},
    {1.0, 0.0, -1.0},   {1.0, 0.0, 1.0},   {0.0, -1.0, -1.0}, {0.0, -1.0, 1.0}, {0.0, 1.0, -1.0},  {0.0, 1.0, 1.0},
    {-1.0, -1.0, -1.0}, {-1.0, -1.0, 1.0}, {-1.0, 1.0, -1.0}, {-1.0, 1.0, 1.0}, {1.0, -1.0, -1.0}, {1.0, -1.0, 1.0},
    {1.0, 1.0, -1.0},   {1.0, 1.0, 1.0},   {0.0, 0.0, 0.0}};

//-- D2Q9 --
template <>
inline constexpr MInt distFld<2, 9>[3] = {4, 4, 0};
template <>
inline constexpr MInt dxQyFld<2, 9> = 3;
template <>
inline constexpr MFloat tp<2, 9>[4] = {F4B9, F1B9, F1B36, 0};

//-- D3Q15 --
template <>
inline constexpr MInt distFld<3, 15>[3] = {6, 8, 0};
template <>
inline constexpr MInt dxQyFld<3, 15> = 5;
template <>
inline constexpr MFloat tp<3, 15>[4] = {F2B9, F1B9, F1B72, 0};

//-- D3Q19 --
template <>
inline constexpr MInt distFld<3, 19>[3] = {6, 12, 0};
template <>
inline constexpr MInt dxQyFld<3, 19> = 5;
template <>
inline constexpr MFloat tp<3, 19>[4] = {F1B3, F1B18, F1B36, 0};

//-- D3Q27 --
template <>
inline constexpr MInt distFld<3, 27>[3] = {6, 12, 8};
template <>
inline constexpr MInt dxQyFld<3, 27> = 9;
template <>
inline constexpr MFloat tp<3, 27>[4] = {F8B27, F2B27, F1B54, F1B216};

} // namespace lbDescriptor


//--wrapper struct for lattice descriptor---------------------------------------
/** \brief  LB lattice descriptor for arrays depending on D
 *  \author Miro Gondrum
 *  \date   27.05.2021
 *  Wrapper struct to easily access LB arrays depending on number of space
 *  dimension D.
 **/
template <MInt D>
struct LbLatticeDescriptorBase {
  LbLatticeDescriptorBase() = delete;

  static constexpr MInt dirFld(MInt i, MInt j, MInt k) {
    if(D == 2) {
      return lbDescriptor::dirFld2[i][j];
    } // else if(D == 3) {
    return lbDescriptor::dirFld3[i][j][k];
  }

  /// Type of the distribution (0:rest, 1:face, 2:edge, 3:corner)
  static constexpr MInt distType(MInt i) {
    MInt sum = 0;
    for(MInt k = 0; k < D; k++) {
      sum += std::abs(idFld(i, k) - 1);
    }
    return sum;
  }
  static constexpr MInt idFld(MInt i, MInt j) { return lbDescriptor::idFld<D>[i][j]; }
  static constexpr MInt intNghbrArray(MInt i, MInt j) { return lbDescriptor::intNghbrArray<D>[i][j]; }
  static constexpr MFloat linearInterpolationCoefficients(MInt i, MInt j) {
    return lbDescriptor::linearInterpolationCoefficients<D>[i][j];
  }
  static constexpr MInt mFld1(MInt i) { return lbDescriptor::mFld1<D>[i]; }
  static constexpr MInt mFld2(MInt i) { return lbDescriptor::mFld2[i]; }

  static constexpr MInt nFld(MInt i, MInt j) { return lbDescriptor::nFld<D>[i][j]; }
  static constexpr MInt pFld(MInt i, MInt j) { return lbDescriptor::pFld<D>[i][j]; }
  static constexpr MInt componentFld(MInt i, MInt j) {
    if(i % 2 == 0) {
      return nFld(i / 2, j);
    } else {
      return pFld(i / 2, j);
    }
  }

  static constexpr MInt nodalConnectivityVector(MInt i, MInt j) {
    return lbDescriptor::nodalConnectivityVector<D>[i][j];
  }
  static constexpr MInt oppositeDist(MInt i) { return lbDescriptor::oppositeDist<D>[i]; }
  static constexpr MFloat ppdfDir(MInt i, MInt j) { return lbDescriptor::ppdfDir<D>[i][j]; }
  static constexpr const MFloat* ppdfDir(MInt i) { return lbDescriptor::ppdfDir<D>[i]; }
};

/** \brief  LB lattice descriptor for arrays depending on D and Q
 *  \author Miro Gondrum
 *  \date   27.05.2021
 *  Wrapper struct to easily access LB arrays depending on number of space
 *  dimension D and number of population Q.
 **/
template <MInt D, MInt Q>
struct LbLatticeDescriptor : public LbLatticeDescriptorBase<D> {
  LbLatticeDescriptor() = delete;

  static constexpr MInt lastId() { return Q - 1; }
  static constexpr MInt d() { return D; }
  static constexpr MInt q() { return Q; }

  static constexpr MInt distFld(MInt i) { return lbDescriptor::distFld<D, Q>[i]; }
  static constexpr MInt dxQyFld() { return lbDescriptor::dxQyFld<D, Q>; }
  static constexpr MFloat tp(MInt i) { return lbDescriptor::tp<D, Q>[i]; }
};

#endif
