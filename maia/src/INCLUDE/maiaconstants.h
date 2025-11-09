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


#ifndef MAIA_CONSTANTS_H
#define MAIA_CONSTANTS_H
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file defines the constants used in MAIA
////////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <limits>
#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"
#include "compiler_config.h"
////////////////////////////////////////////////////////////////////////////////

static constexpr MFloat F0 = 0.0;
static constexpr MFloat F1 = 1.0;
static constexpr MFloat F2 = 2.0;
static constexpr MFloat F3 = 3.0;
static constexpr MFloat F4 = 4.0;
static constexpr MFloat F5 = 5.0;
static constexpr MFloat F6 = 6.0;
static constexpr MFloat F7 = 7.0;
static constexpr MFloat F8 = 8.0;
static constexpr MFloat F9 = 9.0;
static constexpr MFloat F10 = 10.0;
static constexpr MFloat F1B2 = 1.0 / 2.0;
static constexpr MFloat F1B4 = 1.0 / 4.0;
static constexpr MFloat F1B3 = 1.0 / 3.0;
static constexpr MFloat F1B5 = 1.0 / 5.0;
static constexpr MFloat F1B6 = 1.0 / 6.0;
static constexpr MFloat F1B9 = 1.0 / 9.0;
static constexpr MFloat F1B10 = 1.0 / 10.0;
static constexpr MFloat F1B12 = 1.0 / 12.0;
static constexpr MFloat F1B18 = 1.0 / 18.0;
static constexpr MFloat F1B20 = 1.0 / 20.0;
static constexpr MFloat F1B30 = 1.0 / 30.0;
static constexpr MFloat F1B32 = 1.0 / 32.0;
static constexpr MFloat F1B36 = 1.0 / 36.0;
static constexpr MFloat F1B54 = 1.0 / 54.0;
static constexpr MFloat F1B72 = 1.0 / 72.0;
static constexpr MFloat F2B3 = 2.0 / 3.0;
static constexpr MFloat F2B9 = 2.0 / 9.0;
static constexpr MFloat F3B2 = 3.0 / 2.0;
static constexpr MFloat F3B4 = 3.0 / 4.0;
static constexpr MFloat F4B3 = 4.0 / 3.0;
static constexpr MFloat F4B6 = 4.0 / 6.0;
static constexpr MFloat F5B6 = 5.0 / 6.0;
static constexpr MFloat F4B9 = 4.0 / 9.0;
static constexpr MFloat F1BCSsq = 3.0;
static constexpr MFloat F1B8 = 1.0 / 8.0;
static constexpr MFloat F1B80 = 1.0 / 80.0;
static constexpr MFloat F3B8 = 3.0 / 8.0;
static constexpr MFloat F5B8 = 5.0 / 8.0;
static constexpr MFloat F7B3 = 7.0 / 3.0;
static constexpr MFloat F7B8 = 7.0 / 8.0;
static constexpr MFloat F8B3 = 8.0 / 3.0;

static constexpr MFloat F1B64 = 1.0 / 64.0;
static constexpr MFloat F2B27 = 2.0 / 27.0;
static constexpr MFloat F3B64 = 3.0 / 64.0;
static constexpr MFloat F9B8 = 9.0 / 8.0;
static constexpr MFloat F9B64 = 9.0 / 64.0;
static constexpr MFloat F27B64 = 27.0 / 64.0;
static constexpr MFloat F1B19 = 1.0 / 19.0;
static constexpr MFloat F11B2394 = 11.0 / 2394.0;
static constexpr MFloat F1B63 = 1.0 / 63.0;
static constexpr MFloat F4B1197 = 4.0 / 1197.0;
static constexpr MFloat F1B252 = 1.0 / 252.0;
static constexpr MFloat F1B24 = 1.0 / 24.0;
static constexpr MFloat F23B24 = 23.0 / 24.0;
static constexpr MFloat F475B63 = 475.0 / 63.0;
static constexpr MFloat F5B399 = 5.0 / 399.0;
static constexpr MFloat F1B21 = 1.0 / 21.0;
static constexpr MFloat F32B45 = 32.0 / 45.0;
static constexpr MFloat F49B90 = 49.0 / 90.0;

static constexpr MFloat F9B1024 = 9.0 / 1024.0;
static constexpr MFloat F15B1024 = 15.0 / 1024.0;
static constexpr MFloat F25B16 = 25.0 / 16.0;
static constexpr MFloat F25B1024 = 25.0 / 1024.0;
static constexpr MFloat F45B512 = 45.0 / 512.0;
static constexpr MFloat F75B512 = 75.0 / 512.0;
static constexpr MFloat F225B256 = 225.0 / 256.0;
static constexpr MFloat F10000B259 = 10000.0 / 259.0;

static constexpr MFloat F3B16 = 3.0 / 16.0;
static constexpr MFloat F9B16 = 9.0 / 16.0;
static constexpr MFloat F1B16 = 1.0 / 16.0;

static constexpr MFloat F1B216 = 1.0 / 216.0;

static constexpr MFloat F8B27 = 8.0 / 27.0;

static constexpr MFloat PI = 3.14159265358979323844;
static constexpr MFloat PIB2 = 3.14159265358979323844 / 2.0;
static constexpr MFloat LBCS = 0.577350269189625764509; // sqrt(1/3)
static constexpr MFloat F1BCS = 1.73205080756887729352;
static constexpr MFloat CSsq = 1.0 / 3.0;
static constexpr MFloat SQRT2 = 1.41421356237309504880;
static constexpr MFloat SQRT3 = F1BCS;
static constexpr MFloat SQRT5 = 2.236067977499789695409;
static constexpr MFloat SQRT6 = 2.449489742783178098197;
static constexpr MFloat SQRT7 = 2.64575131106459059050;
static constexpr MFloat SQRT[7] = {0.0, 1.0, SQRT2, SQRT3, 2.0, SQRT5, SQRT6};

static constexpr MFloat F1B2mulF1BCSsq = F1B2 * F1BCSsq;
//! A vector withe the values of 2^2 as Int Terms
/*! It is often used so you only have to calulate it once
 */

static constexpr MInt IPOW3[20] = {1,       3,        9,        27,        81,        243,       729,
                                   2187,    6561,     19683,    59049,     177147,    531441,    1594323,
                                   4782969, 14348907, 43046721, 129140163, 387420489, 1162261467};

/* x > 63 is not covered */
constexpr MLong IPOW2(MInt x) {
  ASSERT(x < 64, "IPOW2(x > 63) will give errneous values for 64Bit Longs");
  return 1L << x;
}
constexpr MFloat FPOW2(MInt x) {
  ASSERT(x < 64, "FPOW2(x > 63) will give errneous values for 64Bit Longs");
  return (MFloat)(1L << x);
}
constexpr MFloat FFPOW2(MInt x) {
  ASSERT(x < 64, "FFPOW2(x > 63) will give errneous values for 64Bit Longs");
  return (F1 / (1L << x));
}

/* x > 31 is not covered */
constexpr MLong IPOW4(MInt x) { return 1L << x << x; }
constexpr MFloat FPOW4(MInt x) { return (MFloat)(1L << x << x); }

static constexpr MFloat FPOW3[12] = {1., 3., 9., 27., 81., 243., 729., 2187., 6561., 19683., 59049., 177147.};


static constexpr MFloat FPOW10[11] = {1.,       10.,       100.,       1000.,       10000.,      100000.,
                                      1000000., 10000000., 100000000., 1000000000., 10000000000.};

static constexpr MFloat FTRIG[91] = {
    1.00000000000000000000e+00,
    9.99847695156391269578e-01,
    9.99390827019095762118e-01,
    9.98629534754573833233e-01,
    9.97564050259824197653e-01,
    9.96194698091745545199e-01,
    9.94521895368273289861e-01,
    9.92546151641321983128e-01,
    9.90268068741570361979e-01,
    9.87688340595137770350e-01,
    9.84807753012208020316e-01,
    9.81627183447663975713e-01,
    9.78147600733805688833e-01,
    9.74370064785235245886e-01,
    9.70295726275996472943e-01,
    9.65925826289068312214e-01,
    9.61261695938318894150e-01,
    9.56304755963035435506e-01,
    9.51056516295153531182e-01,
    9.45518575599316846159e-01,
    9.39692620785908427905e-01,
    9.33580426497201742997e-01,
    9.27183854566787424289e-01,
    9.20504853452440374717e-01,
    9.13545457642600866599e-01,
    9.06307787036649936674e-01,
    8.98794046299167037617e-01,
    8.91006524188367898809e-01,
    8.82947592858926988413e-01,
    8.74619707139395741180e-01,
    8.66025403784438707611e-01,
    8.57167300702112333610e-01,
    8.48048096156425956771e-01,
    8.38670567945424050293e-01,
    8.29037572555041624156e-01,
    8.19152044288991798560e-01,
    8.09016994374947451263e-01,
    7.98635510047292829228e-01,
    7.88010753606722014197e-01,
    7.77145961456970901793e-01,
    7.66044443118978013452e-01,
    7.54709580222772014046e-01,
    7.43144825477394244118e-01,
    7.31353701619170570858e-01,
    7.19339800338651191858e-01,
    7.07106781186547572737e-01,
    6.94658370458997365127e-01,
    6.81998360062498476530e-01,
    6.69130606358858237570e-01,
    6.56059028990507275836e-01,
    6.42787609686539362919e-01,
    6.29320391049837501996e-01,
    6.15661475325658291702e-01,
    6.01815023152048267363e-01,
    5.87785252292473137103e-01,
    5.73576436351046048401e-01,
    5.59192903470746793815e-01,
    5.44639035015027195286e-01,
    5.29919264233204900805e-01,
    5.15038074910054377575e-01,
    5.00000000000000111022e-01,
    4.84809620246337114047e-01,
    4.69471562785890861313e-01,
    4.53990499739546859992e-01,
    4.38371146789077459349e-01,
    4.22618261740699441287e-01,
    4.06736643075800208269e-01,
    3.90731128489273937809e-01,
    3.74606593415911959255e-01,
    3.58367949545300379377e-01,
    3.42020143325668823930e-01,
    3.25568154457156755388e-01,
    3.09016994374947451263e-01,
    2.92371704722736769355e-01,
    2.75637355816999163327e-01,
    2.58819045102520739476e-01,
    2.41921895599667896581e-01,
    2.24951054343864920160e-01,
    2.07911690817759453598e-01,
    1.90808995376544915379e-01,
    1.73648177666930414453e-01,
    1.56434465040230924471e-01,
    1.39173100960065687648e-01,
    1.21869343405147489978e-01,
    1.04528463267653456970e-01,
    8.71557427476581381143e-02,
    6.97564737441254550943e-02,
    5.23359562429439664766e-02,
    3.48994967025010802142e-02,
    1.74524064372833763448e-02,
    0.0,
};

static constexpr MFloat faculty[50] = {
    1.00000000000000000000, 1.00000000000000000000, 2.00000000000000000000, 6.00000000000000000000,
    24.0000000000000000000, 120.000000000000000000, 720.000000000000000000, 5040.00000000000000000,
    40320.0000000000000000, 362880.000000000000000, 3628800.00000000000000, 39916800.0000000000000,
    479001600.000000000000, 6227020800.00000000000, 87178291200.0000000000, 1307674368000.00000000,
    20922789888000.0000000, 355687428096000.000000, 6402373705728000.00000, 121645100408832000.000,
    2432902008176640000.00, 51090942171709440000.0, 1.12400072777760768e21, 2.58520167388849766e22,
    6.20448401733239439e23, 1.55112100433309859e25, 4.03291461126605635e26, 1.08888694504183521e28,
    3.04888344611713860e29, 8.84176199373970195e30, 2.65252859812191058e32, 8.22283865417792281e33,
    2.63130836933693530e35, 8.68331761881188649e36, 2.95232799039604140e38, 1.03331479663861449e40,
    3.71993326789901217e41, 1.37637530912263450e43, 5.23022617466601111e44, 2.03978820811974433e46,
    8.15915283247897734e47, 3.34525266131638071e49, 1.40500611775287989e51, 6.04152630633738356e52,
    2.65827157478844876e54, 1.19622220865480194e56, 5.50262215981208894e57, 2.58623241511168180e59,
    1.24139155925360726e61, 6.08281864034267560e62};

// static constexpr MInt ICUBE[16] = {0, 1, 8, 27, 64, 125, 216, 343, 512, 729, 1000, 1331, 1728, 2197, 2744, 3375};

static constexpr MInt hilbertOrders2D[4][4] = {{0, 2, 3, 1}, {0, 1, 3, 2}, {3, 2, 0, 1}, {3, 1, 0, 2}};

static constexpr MInt hilbertConnections2D[4][4] = {{1, 0, 0, 2}, {0, 1, 1, 3}, {3, 2, 2, 0}, {2, 3, 3, 1}};

static constexpr MInt hilbertOrders3D[24][8] = {
    {0, 2, 3, 1, 5, 7, 6, 4}, {0, 4, 5, 1, 3, 7, 6, 2}, {0, 4, 6, 2, 3, 7, 5, 1}, {5, 4, 0, 1, 3, 2, 6, 7},
    {3, 7, 5, 1, 0, 4, 6, 2}, {6, 4, 5, 7, 3, 1, 0, 2}, {0, 2, 6, 4, 5, 7, 3, 1}, {0, 1, 5, 4, 6, 7, 3, 2},
    {6, 2, 0, 4, 5, 1, 3, 7}, {5, 7, 6, 4, 0, 2, 3, 1}, {3, 2, 6, 7, 5, 4, 0, 1}, {3, 2, 0, 1, 5, 4, 6, 7},
    {5, 7, 3, 1, 0, 2, 6, 4}, {6, 7, 3, 2, 0, 1, 5, 4}, {6, 2, 3, 7, 5, 1, 0, 4}, {5, 1, 3, 7, 6, 2, 0, 4},
    {0, 1, 3, 2, 6, 7, 5, 4}, {3, 1, 0, 2, 6, 4, 5, 7}, {6, 4, 0, 2, 3, 1, 5, 7}, {3, 7, 6, 2, 0, 4, 5, 1},
    {5, 4, 6, 7, 3, 2, 0, 1}, {6, 7, 5, 4, 0, 1, 3, 2}, {3, 1, 5, 7, 6, 4, 0, 2}, {5, 1, 0, 4, 6, 2, 3, 7}};

static constexpr MInt hilbertConnections3D[24][8] = {
    {1, 6, 6, 11, 11, 12, 12, 14},    {0, 2, 2, 3, 3, 4, 4, 5},         {16, 1, 1, 18, 18, 19, 19, 20},
    {12, 20, 20, 1, 1, 11, 11, 18},   {11, 19, 19, 12, 12, 1, 1, 21},   {14, 18, 18, 20, 20, 22, 22, 1},
    {7, 0, 0, 8, 8, 9, 9, 10},        {6, 16, 16, 23, 23, 21, 21, 22},  {21, 14, 14, 6, 6, 23, 23, 11},
    {23, 12, 12, 21, 21, 6, 6, 19},   {22, 11, 11, 14, 14, 20, 20, 6},  {4, 10, 10, 0, 0, 3, 3, 8},
    {3, 9, 9, 4, 4, 0, 0, 13},        {18, 21, 21, 19, 19, 16, 16, 12}, {5, 8, 8, 10, 10, 15, 15, 0},
    {20, 23, 23, 22, 22, 14, 14, 16}, {2, 7, 7, 17, 17, 13, 13, 15},    {19, 22, 22, 16, 16, 18, 18, 23},
    {13, 5, 5, 2, 2, 17, 17, 3},      {17, 4, 4, 13, 13, 2, 2, 9},      {15, 3, 3, 5, 5, 10, 10, 2},
    {8, 13, 13, 9, 9, 7, 7, 4},       {10, 17, 17, 15, 15, 5, 5, 7},    {9, 15, 15, 7, 7, 8, 8, 17}};

static constexpr MBool childCode[6][8] = {{0, 1, 0, 1, 0, 1, 0, 1}, {1, 0, 1, 0, 1, 0, 1, 0}, {0, 0, 1, 1, 0, 0, 1, 1},
                                          {1, 1, 0, 0, 1, 1, 0, 0}, {0, 0, 0, 0, 1, 1, 1, 1}, {1, 1, 1, 1, 0, 0, 0, 0}};
static constexpr uint_fast8_t childCodePro[6] = {0b10101010, 0b01010101, 0b11001100,
                                                 0b00110011, 0b11110000, 0b00001111};

static constexpr MInt nghAcrossCell2D[4][8] = {
    {1, 1, 2, 2, 3, 3, 3, 3}, {0, 0, 3, 3, 2, 2, 2, 2}, {3, 3, 0, 0, 1, 1, 1, 1}, {2, 2, 1, 1, 0, 0, 0, 0}};
static constexpr MInt nghAcrossCell3D[8][26] = {
    {1, 1, 2, 2, 4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7},
    {0, 0, 3, 3, 5, 5, 2, 2, 2, 2, 4, 4, 4, 4, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6},
    {3, 3, 0, 0, 6, 6, 1, 1, 1, 1, 7, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5},
    {2, 2, 1, 1, 7, 7, 0, 0, 0, 0, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4},
    {5, 5, 6, 6, 0, 0, 7, 7, 7, 7, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3},
    {4, 4, 7, 7, 1, 1, 6, 6, 6, 6, 0, 0, 0, 0, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2},
    {7, 7, 4, 4, 2, 2, 5, 5, 5, 5, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1},
    {6, 6, 5, 5, 3, 3, 4, 4, 4, 4, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}};

static constexpr MInt nghInside3D[8][6] = {{-1, 1, -1, 2, -1, 4}, {0, -1, -1, 3, -1, 5}, {-1, 3, 0, -1, -1, 6},
                                           {2, -1, 1, -1, -1, 7}, {-1, 5, -1, 6, 0, -1}, {4, -1, -1, 7, 1, -1},
                                           {-1, 7, 4, -1, 2, -1}, {6, -1, 5, -1, 3, -1}};

static constexpr MInt edgeDiagNeighs3D[24][2] = {{0, 2}, {2, 0}, {0, 3}, {3, 0}, {1, 2}, {2, 1}, {1, 3}, {3, 1},
                                                 {4, 2}, {2, 4}, {3, 4}, {4, 3}, {2, 5}, {5, 2}, {3, 5}, {5, 3},
                                                 {0, 4}, {4, 0}, {1, 4}, {4, 1}, {0, 5}, {5, 0}, {1, 5}, {5, 1}};


static constexpr MInt spaceDiagNeighs3D[48][3] = {
    {0, 2, 4}, {0, 4, 2}, {2, 0, 4}, {2, 4, 0}, {4, 0, 2}, {4, 2, 0}, {0, 3, 4}, {0, 4, 3}, {3, 0, 4}, {3, 4, 0},
    {4, 0, 3}, {4, 3, 0}, {1, 2, 4}, {1, 4, 2}, {2, 1, 4}, {2, 4, 1}, {4, 1, 2}, {4, 2, 1}, {1, 3, 4}, {1, 4, 3},
    {3, 1, 4}, {3, 4, 1}, {4, 1, 3}, {4, 3, 1}, {0, 2, 5}, {0, 5, 2}, {2, 0, 5}, {2, 5, 0}, {5, 0, 2}, {5, 2, 0},
    {0, 3, 5}, {0, 5, 3}, {3, 0, 5}, {3, 5, 0}, {5, 0, 3}, {5, 3, 0}, {1, 2, 5}, {1, 5, 2}, {2, 1, 5}, {2, 5, 1},
    {5, 1, 2}, {5, 2, 1}, {1, 3, 5}, {1, 5, 3}, {3, 1, 5}, {3, 5, 1}, {5, 1, 3}, {5, 3, 1}};

static constexpr MInt sharedCornerNeighs2D[4][3] = {{0, 2, 6},  // corner (-1,-1)
                                                    {0, 3, 7},  // corner (-1, 1)
                                                    {1, 2, 8},  // corner (1 ,-1)
                                                    {1, 3, 9}}; // corner (1 , 1)

static constexpr MInt sharedCornerNeighs3D[8][7] = {{0, 2, 4, 6, 10, 14, 18},  // corner (-1,-1,-1)
                                                    {0, 2, 5, 6, 11, 15, 19},  // corner (-1,-1, 1)
                                                    {0, 3, 4, 7, 10, 16, 20},  // corner (-1, 1,-1)
                                                    {0, 3, 5, 7, 11, 17, 21},  // corner (-1, 1, 1)
                                                    {1, 2, 4, 8, 12, 14, 22},  // corner ( 1,-1,-1)
                                                    {1, 2, 5, 8, 13, 15, 23},  // corner ( 1,-1, 1)
                                                    {1, 3, 4, 9, 12, 16, 24},  // corner ( 1, 1,-1)
                                                    {1, 3, 5, 9, 13, 17, 25}}; // corner ( 1, 1, 1)

static constexpr MInt traverseCorners2D[4] = {0, 2, 3, 1};
static constexpr MInt traverseCorners3D[8] = {0, 2, 6, 4, 5, 7, 3, 1};


static constexpr MInt oppositeDirGrid[26] = {1,  0,  3,  2,  5,  4,  9,  8,  7,  6,  13, 12, 11,
                                             10, 17, 16, 15, 14, 25, 24, 23, 22, 21, 20, 19, 18};

static constexpr MInt childPos2D[4] = {6, 5, 7, 4};
static constexpr MInt childPos3D[8] = {18, 22, 20, 24, 19, 23, 21, 25};

//! Max Nr. of Dimensions (3)
static constexpr MInt MAX_SPACE_DIMENSIONS = 3;

//! Epsilon
#ifdef MAIA_PGI_COMPILER
// The PGI compiler uses an old Standard Library, where the epsilon() is not constexpr
static const MFloat MFloatEps = std::numeric_limits<MFloat>::epsilon();
#else
static constexpr MFloat MFloatEps = std::numeric_limits<MFloat>::epsilon();
#endif

//! NaN
static const MFloat MFloatNaN = std::numeric_limits<MFloat>::quiet_NaN();

//! Float max
static const MFloat MFloatMax = std::numeric_limits<MFloat>::max();
//! Sqrt of max floating point value
static const MFloat MFloatMaxSqrt = std::sqrt(std::numeric_limits<MFloat>::max());

// Memory conversion factors
static constexpr MFloat c_byteToGByte = F1 / (1024 * 1024 * 1024);

#endif
