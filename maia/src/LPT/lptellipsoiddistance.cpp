// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lptellipsoiddistance.h"

#include "INCLUDE/maiaconstants.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "UTIL/timer.h"

using namespace std;

/*Subroutine to calculate the distance of closest approach of two arbitrary ellipsoids*/
EllipsoidDistance::EllipsoidDistance(MFloat v1[3], MFloat v2[3], MFloat v3[3], MFloat v4[3], MFloat v5[3], MFloat v6[3],
                                     MFloat v7[3], MFloat s1, MFloat s2, MFloat s3, MFloat s4, MFloat s5, MFloat s6) {
  TRACE();

  for(MInt i = 0; i < 3; i++) {
    di[i] = v1[i];
    l1i[i] = v2[i];
    m1i[i] = v3[i];
    n1i[i] = v4[i];
    l2i[i] = v5[i];
    m2i[i] = v6[i];
    n2i[i] = v7[i];
  }
  a1 = s1;
  b1 = s2;
  c1 = s3;
  a2 = s4;
  b2 = s5;
  c2 = s6;

  pi = 4.0 * atan(1.0);
  gratio = (1.0 + sqrt(5.0)) * 0.5;
  o1g = 1.0 / (1.0 + gratio);
  tolerance = 1.0E-3;
  delt = 1.0E-5;
}

MFloat EllipsoidDistance::ellipsoids(void) {
  TRACE();

  MFloat sol;
  MFloat theta1, theta2, theta3, theta4, f1, f4, f3, f5, f6; // Search algorithm variables

  norm(di, d);
  norm(l1i, l1);
  norm(l2i, l2);
  norm(m1i, m1);
  norm(m2i, m2);
  norm(n1i, n1);
  norm(n2i, n2);

  crossP(d, l1, dxp0);
  if(mag(dxp0) < 1.E-14) crossP(d, m1, dxp0);
  norm(dxp0, p0);
  crossP(p0, d, dxp0);

  theta1 = 0.0;
  theta2 = pi;
  theta3 = theta1 + o1g * (theta2 - theta1);
  theta4 = theta2 - o1g * (theta2 - theta1);
  f1 = plane_int(theta1);
  f3 = plane_int(theta3);
  f4 = plane_int(theta4);

  /////////////////////////////search algorithm loop///////////////////////////////
  while((theta2 - theta1) > tolerance) {
    if((f1 <= f3 && f3 <= f4) || (f3 <= f1 && f1 <= f4)) // case a
    {
      theta1 = theta3;
      f1 = f3;
      theta3 = theta4;
      f3 = f4;
      theta4 = theta2 - o1g * (theta2 - theta1);
      f4 = plane_int(theta4);
    } else if((f1 <= f4 && f4 <= f3) || (f4 <= f1 && f1 <= f3)) // case b
    {
      theta2 = theta4;
      theta4 = theta3;
      f4 = f3;
      theta3 = theta1 + o1g * (theta2 - theta1);
      f3 = plane_int(theta3);
    } else // case c
    {
      f5 = plane_int(theta1 + delt);
      if(f5 > f1) {
        theta2 = theta3;
        theta3 = theta1 + o1g * (theta2 - theta1);
        theta4 = theta2 - o1g * (theta2 - theta1);
        f3 = plane_int(theta3);
        f4 = plane_int(theta4);
      } else {
        f6 = plane_int(theta1 - delt);
        if(f6 < f1) {
          theta2 = theta3;
          theta3 = theta1 + o1g * (theta2 - theta1);
          theta4 = theta2 - o1g * (theta2 - theta1);
          f3 = plane_int(theta3);
          f4 = plane_int(theta4);
        } else {
          theta1 = theta4;
          f1 = f4;
          theta3 = theta1 + o1g * (theta2 - theta1);
          theta4 = theta2 - o1g * (theta2 - theta1);
          f3 = plane_int(theta3);
          f4 = plane_int(theta4);
        }
      }
    }
  }
  sol = f4;


  return sol;
}
////////////////////////////////////////////////////////////////////////////
///////////////////////*FUNCTIONS*///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/*Cross Product of two vectors*/
void EllipsoidDistance::crossP(MFloat x[3], MFloat y[3], MFloat z[3]) {
  TRACE();

  z[0] = x[1] * y[2] - x[2] * y[1];
  z[1] = x[2] * y[0] - x[0] * y[2];
  z[2] = x[0] * y[1] - x[1] * y[0];
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*Normalization of a vector*/
void EllipsoidDistance::norm(MFloat vec[3], MFloat nvec[3]) {
  TRACE();

  MFloat length;
  length = mag(vec);
  if(approx(length, 0.0, MFloatEps)) {
    nvec[0] = 0.0;
    nvec[1] = 0.0;
    nvec[2] = 0.0;
  } else {
    nvec[0] = vec[0] / length;
    nvec[1] = vec[1] / length;
    nvec[2] = vec[2] / length;
  }
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*Magnitude of a vector*/
MFloat EllipsoidDistance::mag(MFloat V[3]) { return sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]); }
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*Dot Product of two vectors*/
MFloat EllipsoidDistance::dotP(MFloat Vec1[3], MFloat Vec2[3]) {
  return Vec1[0] * Vec2[0] + Vec1[1] * Vec2[1] + Vec1[2] * Vec2[2];
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*Ellipses formed with the intersection of the ellipsoid with the plane*/
MFloat EllipsoidDistance::plane_int(MFloat theta) {
  TRACE();

  MFloat alpha1, beta1, gamma1, v1, u1, alpha2, beta2, gamma2, v2, u2;
  MFloat l1x, l1y, m1x, m1y, n1x, n1y, l2x, l2y, m2x, m2y, n2x, n2y;
  MFloat a2d1, a2d2, b2d1, b2d2, angle1, angle2; // Ellipses variables
  MInt j;

  for(j = 0; j < 3; j++)
    p[j] = cos(theta) * p0[j] + sin(theta) * dxp0[j];
  crossP(p, d, s);
  l1x = dotP(l1, d);
  l1y = dotP(l1, s);
  m1x = dotP(m1, d);
  m1y = dotP(m1, s);
  n1x = dotP(n1, d);
  n1y = dotP(n1, s);
  alpha1 = l1x * l1x / (a1 * a1) + m1x * m1x / (b1 * b1) + n1x * n1x / (c1 * c1);
  beta1 = l1y * l1x / (a1 * a1) + m1y * m1x / (b1 * b1) + n1x * n1y / (c1 * c1);
  gamma1 = l1y * l1y / (a1 * a1) + m1y * m1y / (b1 * b1) + n1y * n1y / (c1 * c1);
  v1 = sqrt(4.0 * beta1 * beta1 + (alpha1 - gamma1) * (alpha1 - gamma1));
  angle1 = 0.5 * atan2(-2.0 * beta1, -(alpha1 - gamma1));
  u1 = gamma1 + v1 * sin(angle1) * sin(angle1);
  b2d1 = 1.0 / sqrt(u1);
  a2d1 = 1.0 / sqrt(u1 - v1);
  l2x = dotP(l2, d);
  l2y = dotP(l2, s);
  m2x = dotP(m2, d);
  m2y = dotP(m2, s);
  n2x = dotP(n2, d);
  n2y = dotP(n2, s);
  alpha2 = l2x * l2x / (a2 * a2) + m2x * m2x / (b2 * b2) + n2x * n2x / (c2 * c2);
  beta2 = l2y * l2x / (a2 * a2) + m2y * m2x / (b2 * b2) + n2x * n2y / (c2 * c2);
  gamma2 = l2y * l2y / (a2 * a2) + m2y * m2y / (b2 * b2) + n2y * n2y / (c2 * c2);
  v2 = sqrt(4.0 * beta2 * beta2 + (alpha2 - gamma2) * (alpha2 - gamma2));
  angle2 = 0.5 * atan2(-2.0 * beta2, -(alpha2 - gamma2));
  u2 = gamma2 + v2 * sin(angle2) * sin(angle2);
  b2d2 = 1.0 / sqrt(u2);
  a2d2 = 1.0 / sqrt(u2 - v2);


  return distance2d(a2d1, b2d1, a2d2, b2d2, angle1, angle2);
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*Distance of Closest Approach of two arbitrary ellipses*/
MFloat EllipsoidDistance::distance2d(MFloat a1_input, MFloat b1_input, MFloat a2_input, MFloat b2_input, MFloat angle1,
                                     MFloat angle2) {
  TRACE();

  MFloat eps1, eps2, k1dotd, k2dotd, k1dotk2, nu, Ap[2][2];
  MFloat lambdaplus, lambdaminus, bp2, ap2, cosphi, tanphi2, delta, dp, qu2;
  MFloat A, B, C, D, E, alpha, beta, gamma, P, Q, A2, B2;
  complex<MFloat> U, y, qu, calpha;

  eps1 = sqrt(1.0 - (b1_input * b1_input) / (a1_input * a1_input));
  eps2 = sqrt(1.0 - (b2_input * b2_input) / (a2_input * a2_input));
  k1dotd = cos(angle1);
  k2dotd = cos(angle2);
  k1dotk2 = cos(angle2 - angle1);
  nu = a1_input / b1_input - 1.0;
  Ap[0][0] =
      b1_input * b1_input / (b2_input * b2_input)
      * (1.0 + 0.5 * (1.0 + k1dotk2) * (nu * (2.0 + nu) - eps2 * eps2 * (1.0 + nu * k1dotk2) * (1.0 + nu * k1dotk2)));
  Ap[1][1] =
      b1_input * b1_input / (b2_input * b2_input)
      * (1.0 + 0.5 * (1.0 - k1dotk2) * (nu * (2.0 + nu) - eps2 * eps2 * (1.0 - nu * k1dotk2) * (1.0 - nu * k1dotk2)));
  Ap[0][1] = b1_input * b1_input / (b2_input * b2_input) * 0.5 * sqrt(1.0 - k1dotk2 * k1dotk2)
             * (nu * (2.0 + nu) + eps2 * eps2 * (1.0 - nu * nu * k1dotk2 * k1dotk2));
  lambdaplus =
      0.5 * (Ap[0][0] + Ap[1][1]) + sqrt(0.25 * (Ap[0][0] - Ap[1][1]) * (Ap[0][0] - Ap[1][1]) + Ap[0][1] * Ap[0][1]);
  lambdaminus =
      0.5 * (Ap[0][0] + Ap[1][1]) - sqrt(0.25 * (Ap[0][0] - Ap[1][1]) * (Ap[0][0] - Ap[1][1]) + Ap[0][1] * Ap[0][1]);
  bp2 = 1.0 / sqrt(lambdaplus);
  ap2 = 1.0 / sqrt(lambdaminus);

  if(approx(k1dotk2, 1.0, MFloatEps)) {
    if(Ap[0][0] > Ap[1][1])
      cosphi = b1_input * k1dotd / (a1_input * sqrt(1.0 - eps1 * eps1 * k1dotd * k1dotd));
    else
      cosphi = sqrt(1.0 - k1dotd * k1dotd) / sqrt(1.0 - eps1 * eps1 * k1dotd * k1dotd);
  } else {
    cosphi = 1.0
             / sqrt(2.0 * (Ap[0][1] * Ap[0][1] + (lambdaplus - Ap[0][0]) * (lambdaplus - Ap[0][0]))
                    * (1.0 - eps1 * eps1 * k1dotd * k1dotd))
             * (Ap[0][1] / sqrt(1.0 + k1dotk2)
                    * (b1_input / a1_input * k1dotd + k2dotd + (b1_input / a1_input - 1.0) * k1dotd * k1dotk2)
                + (lambdaplus - Ap[0][0]) / sqrt(1.0 - k1dotk2)
                      * (b1_input / a1_input * k1dotd - k2dotd - (b1_input / a1_input - 1.0) * k1dotd * k1dotk2));
  }

  delta = ap2 * ap2 / (bp2 * bp2) - 1.0;
  if(approx(delta, 0.0, MFloatEps) || approx(cosphi, 0.0, MFloatEps))
    dp = 1.0 + ap2;
  else {
    tanphi2 = 1.0 / (cosphi * cosphi) - 1.0;
    A = -(1.0 + tanphi2) / (bp2 * bp2);
    B = -2.0 * (1.0 + tanphi2 + delta) / bp2;
    C = -tanphi2 - (1.0 + delta) * (1.0 + delta) + (1.0 + (1.0 + delta) * tanphi2) / (bp2 * bp2);
    D = 2.0 * (1.0 + tanphi2) * (1.0 + delta) / bp2;
    E = (1.0 + tanphi2 + delta) * (1.0 + delta);
    A2 = A * A;
    B2 = B * B;
    alpha = -3.0 * B2 / (8.0 * A2) + C / A;
    beta = B2 * B / (8.0 * A2 * A) - B * C / (2.0 * A2) + D / A;
    gamma = -3.0 * B2 * B2 / (256.0 * A2 * A2) + C * B2 / (16.0 * A2 * A) - B * D / (4.0 * A2) + E / A;
    calpha = complex<MFloat>(alpha, 0.0);
    if(approx(beta, 0.0, MFloatEps)) {
      qu = complex<MFloat>(-B / (4.0 * A), 0.0)
           + sqrt(0.5 * (-calpha + sqrt(complex<MFloat>(alpha * alpha - 4.0 * gamma, 0.0))));
    } else {
      P = -F1B12 * alpha * alpha - gamma;
      Q = -alpha * alpha * alpha / 108.0 + F1B3 * alpha * gamma - 0.125 * beta * beta;
      U = c_cbrt(complex<MFloat>(-Q * 0.5, 0.0) + sqrt(complex<MFloat>(Q * Q * 0.25 + P * P * P / 27.0, 0.0)));
      if(approx(real(U), 0.0, MFloatEps) && approx(imag(U), 0.0, MFloatEps))
        y = -5.0 * F1B6 * calpha - c_cbrt(complex<MFloat>(Q, 0.0));
      else
        y = -5.0 * F1B6 * calpha + U - complex<MFloat>(P, 0.0) / (3.0 * U);
      qu = complex<MFloat>(-B / (4.0 * A), 0.0)
           + 0.5
                 * (sqrt(calpha + 2.0 * y)
                    + sqrt(-(3.0 * calpha + 2.0 * y + 2.0 * complex<MFloat>(beta, 0.0) / sqrt(calpha + 2.0 * y))));
    }
    qu2 = abs(qu * qu);
    dp = sqrt((qu2 - 1.0) / delta * (1.0 + bp2 * (1.0 + delta)) * (1.0 + bp2 * (1.0 + delta)) / qu2
              + (1.0 - (qu2 - 1.0) / delta) * (1.0 + bp2) * (1.0 + bp2) / qu2);
  }


  return dp * b1_input / sqrt(1.0 - eps1 * eps1 * k1dotd * k1dotd);
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/* Principal cubic root of a complex number */
complex<MFloat> EllipsoidDistance::c_cbrt(complex<MFloat> x) {
  TRACE();

  MFloat a, b, r, phi, rn;
  a = real(x);
  b = imag(x);
  r = sqrt(a * a + b * b);
  phi = atan2(b, a);
  phi /= 3.0;
  rn = cbrt(r);

  return complex<MFloat>(rn * cos(phi), rn * sin(phi));
}
