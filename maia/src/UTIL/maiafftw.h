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


#ifndef MAIA_FFTW_H
#define MAIA_FFTW_H

#include <complex>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <random>

#include "../COMM/mpioverride.h"
#include "../INCLUDE/maiaconstants.h"
#include "../INCLUDE/maiatypes.h"
#include "../IO/context.h"
#include "../MEMORY/scratch.h"
#include "debug.h"

namespace maia {
namespace math {

inline MInt getGlobalPosFFTW(MInt i0, MInt i1, MInt i2, MInt ny, MInt nz) { return i2 + nz * (i1 + ny * i0); }

inline MFloat triadImag(fftw_complex& a, fftw_complex& b, fftw_complex& c) {
  return (a[0] * b[1] * c[0] + a[1] * b[0] * c[0] - a[0] * b[0] * c[1] + a[1] * b[1] * c[1]);
}


/** \brief Generates a velocity field from Fourier-modes using FFTW
 *
 * \author lennart, corrected version of george's version in lb solver
 *
 * discrete spectrum is computed on unity cube, no spatial scaling required
 * velocity field is computed for u_rms = 1, hence uPhysField subsequently has to be scaled by the magnitude of the
 * fluctuations, e.g. uPhysField *= m_UInfinity
 *
 * (u,v,w)PhysField: complex velocity in physical space
 *
 * \param uPhysField pointer to a fftw_complex for u
 * \param vPhysField pointer to a fftw_complex for v
 * \param wPhysField pointer to a fftw_complex for w
 * \param kpRatio ratio of peak wave number to minimum wave number, kp/k0
 * \param spectrumId prescribed energy spectrum
 */

inline MUlong initFft(fftw_complex*& uPhysField, fftw_complex*& nabla2P, MInt nx, MInt ny, MInt nz,
                      const MFloat kpRatio, const MInt spectrumId, MIntScratchSpace& fftInfo, const MPI_Comm comm,
                      const MBool computeNabla2P = false) {
  DEBUG("initFft entry", MAIA_DEBUG_TRACE_IN);

  const MInt size = nx * ny * nz;
  const ptrdiff_t rank = 3;
  const MInt howmany = 3;
  if(nx % 2 != 0 || ny % 2 != 0 || nz % 2 != 0) {
    std::stringstream errorMessage;
    errorMessage << " FFTInit: domainsize must NOT be an odd number! " << nx << " " << ny << " " << nz << std::endl;
    mTerm(1, AT_, errorMessage.str());
  }

  m_log << " --- initializing FFTW --- " << std::endl;
  m_log << " domain size = " << nx << "x" << ny << "x" << nz << std::endl;

  MInt domainId;
  MInt noDomains;
  MPI_Comm_rank(comm, &domainId);
  MPI_Comm_size(comm, &noDomains);
  MPI_Comm MPI_COMM_FFTW;
  MInt maxRank = 1;
  maxRank = mMin(nz, noDomains);
  while(nz % maxRank != 0) {
    maxRank--;
  }

  MInt color = (domainId < maxRank) ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(comm, color, domainId, &MPI_COMM_FFTW, AT_, "MPI_COMM_FFTW");
  m_log << "no ranks for fft: " << maxRank << std::endl;

  MInt tmpHowMany = 2;
  ptrdiff_t alloc_local, local_n0, local_0_start;
  ptrdiff_t alloc_localP, local_n0P, local_0_startP;
  ptrdiff_t alloc_localTmp, tmpLocal_n0, tmpLocal_0_start;

  const ptrdiff_t n[3] = {nx, ny, nz};
#ifndef MAIA_WINDOWS
  alloc_local = (domainId < maxRank) ? fftw_mpi_local_size_many(rank, n, howmany, FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_FFTW,
                                                                &local_n0, &local_0_start)
                                     : 1;
  alloc_localTmp = (domainId < maxRank) ? fftw_mpi_local_size_many(rank, n, tmpHowMany, FFTW_MPI_DEFAULT_BLOCK,
                                                                   MPI_COMM_FFTW, &tmpLocal_n0, &tmpLocal_0_start)
                                        : 1;
  alloc_localP = (domainId < maxRank) ? fftw_mpi_local_size_many(rank, n, 1, FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_FFTW,
                                                                 &local_n0P, &local_0_startP)
                                      : 1;
#endif

  fftInfo[0] = maxRank;
  fftInfo[1] = (domainId < maxRank) ? local_n0 : 0;
  fftInfo[2] = (domainId < maxRank) ? local_0_start : size;
  fftInfo[3] = (domainId < maxRank) ? alloc_local : 3;

  /*! \page propertyPage1
    \section  referenceCubeSize
    <code>MInt MAIAMath::referenceCubeSize</code>\n
    default = <code>1</code>\n\n
    Defines the length of the box in which the FFT is used. \n
    On this account, the minimal wavenumber is defined by kmin=2*pi/referenceCubeSize\n
    Possible values are:
    <ul>
    <li>positive integers </li>
    </ul>
    Keywords: <i>FFT, MATH</i>
  */
  static const MInt referenceCubeSize =
      ((Context::propertyExists("referenceCubeSize")) ? Context::getBasicProperty<MInt>("referenceCubeSize", AT_) : 1);

  static const MInt kMinSpec = ((Context::propertyExists("kMin")) ? Context::getBasicProperty<MInt>("kMin", AT_) : 0);

  static const MInt kMaxSpec = ((Context::propertyExists("kMax")) ? Context::getBasicProperty<MInt>("kMax", AT_) : 0);

  static const int64_t seed0 = (int64_t)(
      (Context::propertyExists("randomDeviceSeed")) ? Context::getBasicProperty<MFloat>("randomDeviceSeed", AT_) : -1);

  static MUlong seed = (unsigned)((seed0 > -1) ? seed0 : std::random_device()());
  MPI_Bcast(&seed, 1, MPI_UNSIGNED_LONG, 0, comm, AT_, "seed");
  std::mt19937_64 gen(seed);
  std::normal_distribution<> distr(F0, F1);

  uPhysField = (fftw_complex*)fftw_malloc(alloc_local * sizeof(fftw_complex));

  if(computeNabla2P) nabla2P = (fftw_complex*)fftw_malloc(alloc_localP * sizeof(fftw_complex));

  if(domainId < maxRank) {
    ptrdiff_t rSlice_0_start = nx - (MInt)local_0_start - (MInt)local_n0 + 1;
    ScratchSpace<fftw_complex> uHatField(alloc_local, AT_, "uHatField");
    ScratchSpace<fftw_complex> vHatField(alloc_local, AT_, "vHatField");
    ScratchSpace<fftw_complex> wHatField(alloc_local, AT_, "wHatField");
    ScratchSpace<fftw_complex> ruHatField(alloc_local, AT_, "ruHatField");
    ScratchSpace<fftw_complex> rvHatField(alloc_local, AT_, "rvHatField");
    ScratchSpace<fftw_complex> rwHatField(alloc_local, AT_, "rwHatField");
    ScratchSpace<fftw_complex> tmpHatFieldMem(alloc_localTmp, AT_, "tmpHatFieldMem");

    fftw_complex* hatField = uPhysField;
    fftw_complex* tmpHat = &tmpHatFieldMem[0];

#ifndef MAIA_WINDOWS
    fftw_plan plan = fftw_mpi_plan_many_dft(rank, n, howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, hatField,
                                            hatField, MPI_COMM_FFTW, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan tmpPlan = fftw_mpi_plan_many_dft(rank, n, tmpHowMany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                               tmpHat, tmpHat, MPI_COMM_FFTW, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
    fftw_plan plan;
    fftw_plan tmpPlan;
#endif

    m_log << "random device seed is " << seed << std::endl;
    const MInt refSize =
        (nx <= 256) ? 256 : nx; // This won´t work if you want to compare LES to grids with more than 256^3 cells
    for(MInt k0 = -refSize / 2 + 1; k0 <= refSize / 2; k0++) {
      for(MInt k1 = -refSize / 2 + 1; k1 <= refSize / 2; k1++) {
        for(MInt k2 = -refSize / 2 + 1; k2 <= refSize / 2; k2++) {
          if(k0 < -nx / 2 + 1 || k0 > nx / 2 || k1 < -ny / 2 + 1 || k1 > ny / 2 || k2 < -nz / 2 + 1 || k2 > nz / 2) {
            (void)distr(gen);
            (void)distr(gen);
            (void)distr(gen);
            (void)distr(gen);
            (void)distr(gen);
            (void)distr(gen);
            if(nx == refSize) {
              mTerm(1, AT_, "wrong refSize.");
            }
            continue;
          }

          MInt i0 = (k0 >= 0) ? k0 : k0 + nx;
          MInt i1 = (k1 >= 0) ? k1 : k1 + ny;
          MInt i2 = (k2 >= 0) ? k2 : k2 + nz;
          MInt pos = getGlobalPosFFTW(i0, i1, i2, ny, nz);
          MInt localrPos = pos - (((MInt)rSlice_0_start) * ny * nz);
          MInt localPos = pos - (((MInt)local_0_start) * ny * nz);

          // Collect random numbers on the current FFTW slice
          if(pos >= ((MInt)local_0_start) * nz * ny && pos < ((MInt)(local_0_start + local_n0) * ny * nz)) {
            if(localPos >= ny * nz * local_n0 || localPos < 0) {
              mTerm(1, AT_, "index exceeds array(1a)");
            }
            uHatField[localPos][0] = distr(gen);
            uHatField[localPos][1] = distr(gen);
            vHatField[localPos][0] = distr(gen);
            vHatField[localPos][1] = distr(gen);
            wHatField[localPos][0] = distr(gen);
            wHatField[localPos][1] = distr(gen);
            continue;
          }
          if(pos >= ((MInt)rSlice_0_start) * nz * ny && pos < ((MInt)(rSlice_0_start + local_n0) * ny * nz) && i0 != 0
             && i0 != nx / 2) {
            localrPos = pos - (((MInt)rSlice_0_start) * ny * nz);
            if(localrPos >= ny * nz * local_n0 || localrPos < 0) {
              mTerm(1, AT_, "index exceeds array(1c)");
            }
            ruHatField[localrPos][0] = distr(gen);
            ruHatField[localrPos][1] = distr(gen);
            rvHatField[localrPos][0] = distr(gen);
            rvHatField[localrPos][1] = distr(gen);
            rwHatField[localrPos][0] = distr(gen);
            rwHatField[localrPos][1] = distr(gen);
            continue;
          }
          if(!(pos >= ((MInt)local_0_start) * nz * ny && pos < ((MInt)(local_0_start + local_n0) * ny * nz))
             && !(pos >= ((MInt)rSlice_0_start) * nz * ny && pos < ((MInt)(rSlice_0_start + local_n0) * ny * nz)
                  && i0 != 0 && i0 != nx / 2)) {
            for(MInt k = 0; k < howmany; k++) {
              (void)distr(gen);
              (void)distr(gen);
            }
          }
        }
      }
    }

    const MFloat kMin = F2 * PI / referenceCubeSize;
    const MFloat kp = kpRatio * kMin * referenceCubeSize;

    MFloat eng = F0;
    MFloat eps = F0;
    MFloat length = F0;

    for(MInt i0 = (MInt)local_0_start; i0 < (MInt)(local_0_start + local_n0); i0++) {
      for(MInt i1 = 0; i1 < ny; i1++) {
        for(MInt i2 = 0; i2 < nz; i2++) {
          MFloat k[3];
          MInt j0 = (i0 > nx / 2) ? i0 - nx : i0;
          MInt j1 = (i1 > ny / 2) ? i1 - ny : i1;
          MInt j2 = (i2 > nz / 2) ? i2 - nz : i2;

          /*To limit the spectrum to a given number of wavenumbers.
          kMinSpec and kMaxSpec are integers. Tested with spectrumId = 3*/
          if(kMinSpec > 0 && kMaxSpec > 0) {
            if(j0 > kMaxSpec || j0 < kMinSpec) j0 = 0;
            if(j1 > kMaxSpec || j1 < kMinSpec) j1 = 0;
            if(j2 > kMaxSpec || j2 < kMinSpec) j2 = 0;
          }
          k[0] = ((MFloat)j0) * kMin;
          k[1] = ((MFloat)j1) * kMin;
          k[2] = ((MFloat)j2) * kMin;

          MInt pos = getGlobalPosFFTW(i0, i1, i2, ny, nz);
          MInt rpos = ((nz - i2) % nz) + nz * (((ny - i1) % ny) + ny * ((nx - i0) % nx));
          MInt localPos = pos - (((MInt)local_0_start) * ny * nz);
          MInt localrPos = rpos - (((MInt)rSlice_0_start) * ny * nz);
          ASSERT(pos > -1 && pos < size && rpos > -1 && rpos < size, "");

          MFloat kAbs = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

          MFloat energy = F0;

          if(kAbs < 1e-15) {
            hatField[howmany * localPos][0] = F0;
            hatField[howmany * localPos][1] = F0;
            hatField[howmany * localPos + 1][0] = F0;
            hatField[howmany * localPos + 1][1] = F0;
            hatField[howmany * localPos + 2][0] = F0;
            hatField[howmany * localPos + 2][1] = F0;
            continue;
          }

          // set spectral distribution
          switch(spectrumId) {
            case 0: {
              energy = pow(kAbs / kp, 4.0) * exp(-2.0 * (kAbs / kp) * (kAbs / kp));
              break;
            }
            case 1: {
              energy = POW2(kAbs / kp) * exp(-F1B2 * POW2(kAbs / kp));
              break;
            }
            case 2: {
              // see Schumann and Patterson
              energy = F3B2 * (kAbs / POW2(kp)) * exp(-kAbs / kp);
              break;
            }
            default: {
              mTerm(1, AT_, "Unknown spectrum.");
              break;
            }
          }

          if(j0 > 0 && j1 == 0 && j2 == 0) {
            const MFloat dk = kMin;
            eng += dk * energy;
            eps += F2 * dk * POW2(k[0]) * energy;
            length += F1B2 * PI * dk * energy / k[0];
          }

          // assemble uHat; see Orszag, Phys. Fluids (1969)
          MFloat r[3], s[3];
          const MFloat fac = PI * sqrt(energy) / (sqrt(referenceCubeSize * F2) * referenceCubeSize * kAbs);

          if(i0 == 0 || i0 == nx / 2) {
            localrPos = rpos - (((MInt)local_0_start) * ny * nz);
            r[0] = fac * (uHatField[localPos][0] + uHatField[localrPos][0]);
            s[0] = fac * (uHatField[localPos][1] - uHatField[localrPos][1]);
            r[1] = fac * (vHatField[localPos][0] + vHatField[localrPos][0]);
            s[1] = fac * (vHatField[localPos][1] - vHatField[localrPos][1]);
            r[2] = fac * (wHatField[localPos][0] + wHatField[localrPos][0]);
            s[2] = fac * (wHatField[localPos][1] - wHatField[localrPos][1]);
          } else {
            if(localPos >= ny * nz * local_n0 || localrPos >= ny * nz * local_n0 || localPos < 0 || localrPos < 0) {
              mTerm(1, AT_, "index exceeds array(2)");
            }
            r[0] = fac * (uHatField[localPos][0] + ruHatField[localrPos][0]);
            s[0] = fac * (uHatField[localPos][1] - ruHatField[localrPos][1]);
            r[1] = fac * (vHatField[localPos][0] + rvHatField[localrPos][0]);
            s[1] = fac * (vHatField[localPos][1] - rvHatField[localrPos][1]);
            r[2] = fac * (wHatField[localPos][0] + rwHatField[localrPos][0]);
            s[2] = fac * (wHatField[localPos][1] - rwHatField[localrPos][1]);
          }


          hatField[howmany * localPos][0] = (1.0 - k[0] * k[0] / (kAbs * kAbs)) * r[0]
                                            - k[0] * k[1] / (kAbs * kAbs) * r[1] - k[0] * k[2] / (kAbs * kAbs) * r[2];
          hatField[howmany * localPos][1] = (1.0 - k[0] * k[0] / (kAbs * kAbs)) * s[0]
                                            - k[0] * k[1] / (kAbs * kAbs) * s[1] - k[0] * k[2] / (kAbs * kAbs) * s[2];

          hatField[howmany * localPos + 1][0] = -k[1] * k[0] / (kAbs * kAbs) * r[0]
                                                + (1.0 - k[1] * k[1] / (kAbs * kAbs)) * r[1]
                                                - k[1] * k[2] / (kAbs * kAbs) * r[2];
          hatField[howmany * localPos + 1][1] = -k[1] * k[0] / (kAbs * kAbs) * s[0]
                                                + (1.0 - k[1] * k[1] / (kAbs * kAbs)) * s[1]
                                                - k[1] * k[2] / (kAbs * kAbs) * s[2];

          hatField[howmany * localPos + 2][0] = -k[2] * k[0] / (kAbs * kAbs) * r[0] - k[2] * k[1] / (kAbs * kAbs) * r[1]
                                                + (1.0 - k[2] * k[2] / (kAbs * kAbs)) * r[2];
          hatField[howmany * localPos + 2][1] = -k[2] * k[0] / (kAbs * kAbs) * s[0] - k[2] * k[1] / (kAbs * kAbs) * s[1]
                                                + (1.0 - k[2] * k[2] / (kAbs * kAbs)) * s[2];
        }
      }
    }
    if(computeNabla2P) {
      for(MInt i = 0; i < local_n0 * nz * ny; i++) {
        if(i >= alloc_local / 3) mTerm(1, AT_, "index exceeds array(6)");
        nabla2P[i][0] = F0;
        nabla2P[i][1] = F0;
      }


      for(MInt i = 0; i < 3; i++) {
        for(MInt j = i; j < 3; j++) {
          for(MInt i0 = (MInt)tmpLocal_0_start; i0 < (MInt)(tmpLocal_0_start + tmpLocal_n0); i0++) {
            for(MInt i1 = 0; i1 < ny; i1++) {
              for(MInt i2 = 0; i2 < nz; i2++) {
                MFloat k[3];
                MInt j0 = (i0 > nx / 2) ? i0 - nx : i0;
                MInt j1 = (i1 > ny / 2) ? i1 - ny : i1;
                MInt j2 = (i2 > nz / 2) ? i2 - nz : i2;
                k[0] = ((MFloat)j0) * kMin;
                k[1] = ((MFloat)j1) * kMin;
                k[2] = ((MFloat)j2) * kMin;
                MInt pos = getGlobalPosFFTW(i0, i1, i2, ny, nz);
                MInt localPos = pos - (((MInt)tmpLocal_0_start) * ny * nz);
                if(localPos * tmpHowMany + tmpHowMany - 1 >= alloc_localTmp) mTerm(1, AT_, "index exceeds array(3)");
                const MInt dimIndex0 = (i == 0) ? 0 : ((i == 1) ? 1 : 2);
                const MInt dimIndex1 = (j == 0) ? 0 : ((j == 1) ? 1 : 2);
                tmpHat[localPos * tmpHowMany][0] = -k[j] * hatField[howmany * localPos + dimIndex0][1];
                tmpHat[localPos * tmpHowMany][1] = k[j] * hatField[howmany * localPos + dimIndex0][0];
                tmpHat[localPos * tmpHowMany + 1][0] = -k[i] * hatField[howmany * localPos + dimIndex1][1];
                tmpHat[localPos * tmpHowMany + 1][1] = k[i] * hatField[howmany * localPos + dimIndex1][0];
              }
            }
          }

          fftw_execute(tmpPlan);
          MFloat fac = (i == j) ? F1 : F2;

          for(MInt i0 = (MInt)tmpLocal_0_start; i0 < (MInt)(tmpLocal_0_start + tmpLocal_n0); i0++) {
            for(MInt i1 = 0; i1 < ny; i1++) {
              for(MInt i2 = 0; i2 < nz; i2++) {
                MInt pos = getGlobalPosFFTW(i0, i1, i2, ny, nz);
                MInt localPos = pos - (((MInt)tmpLocal_0_start) * ny * nz);
                if(localPos * tmpHowMany + tmpHowMany - 1 >= alloc_localTmp) mTerm(1, AT_, "index exceeds array(4)");
                if(localPos >= alloc_local / 3) mTerm(1, AT_, "index exceeds array(5)");
                nabla2P[localPos][0] += fac * tmpHat[localPos * tmpHowMany][0] * tmpHat[localPos * tmpHowMany + 1][0];
                tmpHat[localPos * tmpHowMany][0] = 0;
                tmpHat[localPos * tmpHowMany][1] = 0;
                tmpHat[localPos * tmpHowMany + 1][0] = 0;
                tmpHat[localPos * tmpHowMany + 1][1] = 0;
              }
            }
          }
        }
      }

      fftw_destroy_plan(tmpPlan);

      for(MInt i = 0; i < local_n0 * nz * ny; i++) {
        if(i >= alloc_local / 3) mTerm(1, AT_, "index exceeds array(6)");
        nabla2P[i][1] = F0;
      }

      fftw_plan planP = fftw_mpi_plan_many_dft(rank, n, 1, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, nabla2P,
                                               nabla2P, MPI_COMM_FFTW, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(planP);

      for(MInt k = 0; k < local_n0 * nz * ny; k++) {
        if(k >= alloc_local / 3) mTerm(1, AT_, "index exceeds array(7)");
        nabla2P[k][0] /= ((MFloat)size);
        nabla2P[k][1] /= ((MFloat)size);
      }

      fftw_destroy_plan(planP);
      for(MInt i0 = (MInt)local_0_start; i0 < (MInt)(local_0_start + local_n0); i0++) {
        for(MInt i1 = 0; i1 < ny; i1++) {
          for(MInt i2 = 0; i2 < nz; i2++) {
            MFloat k[3];
            MInt j0 = (i0 > nx / 2) ? i0 - nx : i0;
            MInt j1 = (i1 > ny / 2) ? i1 - ny : i1;
            MInt j2 = (i2 > nz / 2) ? i2 - nz : i2;
            k[0] = ((MFloat)j0) * kMin;
            k[1] = ((MFloat)j1) * kMin;
            k[2] = ((MFloat)j2) * kMin;
            MInt pos = getGlobalPosFFTW(i0, i1, i2, ny, nz);
            MInt localPos = pos - (((MInt)local_0_start) * ny * nz);
            MFloat k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            if(localPos >= alloc_local / 3) mTerm(1, AT_, "index exceeds array(8)");
            if(k2 < 1e-15) {
              nabla2P[localPos][0] = F0;
              nabla2P[localPos][1] = F0;
              continue;
            }
            nabla2P[localPos][0] = nabla2P[localPos][0] / k2;
            nabla2P[localPos][1] = nabla2P[localPos][1] / k2;
          }
        }
      }

#ifndef MAIA_WINDOWS
      fftw_plan planP2 = fftw_mpi_plan_many_dft(rank, n, 1, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, nabla2P,
                                                nabla2P, MPI_COMM_FFTW, FFTW_BACKWARD, FFTW_ESTIMATE);
      fftw_execute(planP2);
      fftw_destroy_plan(planP2);
#endif
    }
    MPI_Allreduce(MPI_IN_PLACE, &eng, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE", "eng");
    MPI_Allreduce(MPI_IN_PLACE, &eps, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE", "eps");
    MPI_Allreduce(MPI_IN_PLACE, &length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE", "length");

    if(domainId == 0)
      std::cerr << "SPECTRUM: kinetic energy: " << eng << " (theoretical: 1.5)"
                << ", dissipation rate: " << eps << " (theoretical: " << 72.0 * POW2(kpRatio * PI) << ")"
                << ", integral length/unit length: " << length << " (theoretical: " << F3B8 / kpRatio << ")"
                << std::endl;

    fftw_execute(plan);
  }
  m_log << "fft finished" << std::endl;

  return seed;

  DEBUG("initFft return", MAIA_DEBUG_TRACE_OUT);
}

/* brief returns a normal distributed random-
 * number with mu=mean and sigma=standard deviation
 *
 */
inline MFloat randnormal(MFloat mu, MFloat sigma) {
  TRACE();

  static MBool deviateAvailable = false; //        flag
  static float storedDeviate;            //        deviate from previous calculation
  MFloat polar, rsquared, var1, var2;

  //        If no deviate has been stored, the polar Box-Muller transformation is
  //        performed, producing two independent normally-distributed random
  //        deviates.  One is stored for the next round, and one is returned.
  if(!deviateAvailable) {
    //        choose pairs of uniformly distributed deviates, discarding those
    //        that don't fall within the unit circle
    do {
      var1 = 2.0 * (MFloat(rand()) / MFloat(RAND_MAX)) - 1.0;
      var2 = 2.0 * (MFloat(rand()) / MFloat(RAND_MAX)) - 1.0;
      rsquared = var1 * var1 + var2 * var2;
    } while(rsquared >= 1.0 || approx(rsquared, 0.0, MFloatEps));

    //        calculate polar tranformation for each deviate
    polar = sqrt(-2.0 * log(rsquared) / rsquared);

    //        store first deviate and set flag
    storedDeviate = var1 * polar;
    deviateAvailable = true;

    //        return second deviate

    return var2 * polar * sigma + mu;
  }

  //        If a deviate is available from a previous call to this function, it is
  //        returned, and the flag is set to false.
  else {
    deviateAvailable = false;

    return storedDeviate * sigma + mu;
  }
}

/** brief Generates a single std::complex coefficient of Fourier series
 *
 *  for a given wavenumber k, and a certain energy spectrum
 *  (see Appendix of Orszag, 1969)
 *
 */
inline std::complex<MFloat>* getFourierCoefficient(MFloat* k, MFloat k0, const MFloat Ma) {
  TRACE();

  MFloat r[6], s[6], kAbs, energy;
  std::complex<MFloat> uHat, vHat, wHat;
  std::complex<MFloat>* fourierCoefficient;

  kAbs = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  // the zero-frequency component is always set to zero, so there is no offset
  if(approx(kAbs, 0.0, MFloatEps)) {
    fourierCoefficient = new std::complex<MFloat>[3];

    fourierCoefficient[0] = std::complex<MFloat>(0, 0);
    fourierCoefficient[1] = std::complex<MFloat>(0, 0);
    fourierCoefficient[2] = std::complex<MFloat>(0, 0);


    return fourierCoefficient;

  } else {
    // energy = (kAbs/k0)*(kAbs/k0)*(kAbs/k0)*(kAbs/k0) * exp(-2.0*(kAbs/k0)*(kAbs/k0));
    // energy = pow(kAbs/k0,8.0) * exp(-4.0*(kAbs/k0)*(kAbs/k0)); // set spectral distribution
    energy = pow(kAbs / k0, 4.0) * exp(-2.0 * (kAbs / k0) * (kAbs / k0)); // set spectral distribution
    energy *= exp(2.0) * 0.499 * (Ma * LBCS)
              * (Ma * LBCS); // set maximal fluctuation amplitude to 20% of the freestream velocity (for 128^3: 0.88)

    // determine Fourier coefficients:
    // r and s are Independant random vector fields with independant
    // components (zero mean and rms according to energy spectrum).
    // Each vector has three components for k and another three for -k.

    for(MInt i = 0; i < 6; i++) {
      r[i] = randnormal(0.0, PI * sqrt(energy) / (SQRT2 * kAbs));
      // r[i] = randNumGen.randNorm(0.0, PI*sqrt(energy)/(SQRT2*kAbs));
      s[i] = randnormal(0.0, PI * sqrt(energy) / (SQRT2 * kAbs));
      // s[i] = randNumGen.randNorm(0.0, PI*sqrt(energy)/(SQRT2*kAbs));
    }

    uHat = (1.0 - k[0] * k[0] / (kAbs * kAbs)) * std::complex<MFloat>(r[0] + r[3], s[0] - s[3])
           - k[0] * k[1] / (kAbs * kAbs) * std::complex<MFloat>(r[1] + r[4], s[1] - s[4])
           - k[0] * k[2] / (kAbs * kAbs) * std::complex<MFloat>(r[2] + r[5], s[2] - s[5]);


    vHat = -k[1] * k[0] / (kAbs * kAbs) * std::complex<MFloat>(r[0] + r[3], s[0] - s[3])
           + (1.0 - k[1] * k[1] / (kAbs * kAbs)) * std::complex<MFloat>(r[1] + r[4], s[1] - s[4])
           - k[1] * k[2] / (kAbs * kAbs) * std::complex<MFloat>(r[2] + r[5], s[2] - s[5]);


    wHat = -k[2] * k[0] / (kAbs * kAbs) * std::complex<MFloat>(r[0] + r[3], s[0] - s[3])
           - k[2] * k[1] / (kAbs * kAbs) * std::complex<MFloat>(r[1] + r[4], s[1] - s[4])
           + (1.0 - k[2] * k[2] / (kAbs * kAbs)) * std::complex<MFloat>(r[2] + r[5], s[2] - s[5]);

    fourierCoefficient = new std::complex<MFloat>[3];

    // fourierCoefficient[0] = std::complex<MFloat>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);// uHat;
    // fourierCoefficient[1] = std::complex<MFloat>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);//vHat;
    // fourierCoefficient[2] = std::complex<MFloat>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);//wHat;

    fourierCoefficient[0] = uHat;
    fourierCoefficient[1] = vHat;
    fourierCoefficient[2] = wHat;
    return fourierCoefficient;
  }
}


/** brief Generates a single std::complex coefficient of Fourier series
 *
 *  for a given wavenumber k, and a certain energy spectrum
 *
 */
inline std::complex<MFloat>* getFourierCoefficient2(MFloat* k, MFloat k0) {
  static std::mt19937 randNumGen;
  static std::uniform_real_distribution<> distrib{0.0, 1.0};

  MFloat r[3], kAbs, energy;
  std::complex<MFloat> uHat, vHat, wHat;
  std::complex<MFloat>* fourierCoefficient;

  kAbs = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  // the zero-frequency component is always set to zero, so there is no offset
  if(approx(kAbs, 0.0, MFloatEps)) {
    fourierCoefficient = new std::complex<MFloat>[3];

    fourierCoefficient[0] = std::complex<MFloat>(0, 0);
    fourierCoefficient[1] = std::complex<MFloat>(0, 0);
    fourierCoefficient[2] = std::complex<MFloat>(0, 0);


    return fourierCoefficient;

  } else {
    // energy = pow(kAbs/k0,4.0) * exp(-2.0*(kAbs/k0)*(kAbs/k0));
    // energy = ( m_Ma*LBCS * m_Ma*LBCS ) * pow(kAbs/k0,4.0) * exp(-2.0*(kAbs/k0)*(kAbs/k0));
    // energy = 0.135335;
    energy = (kAbs / k0) * (kAbs / k0) * (kAbs / k0) * (kAbs / k0) * exp(-2.0 * (kAbs / k0) * (kAbs / k0));

    // determine Fourier coefficients:
    for(MInt i = 0; i < 3; i++) {
      r[i] = 2.0 * PI * distrib(randNumGen);
    }

    uHat = std::complex<MFloat>(cos(r[0]) * energy, sin(r[0]) * energy);
    vHat = std::complex<MFloat>(cos(r[1]) * energy, sin(r[1]) * energy);
    wHat = std::complex<MFloat>(cos(r[2]) * energy, sin(r[2]) * energy);

    fourierCoefficient = new std::complex<MFloat>[3];

    fourierCoefficient[0] = uHat;
    fourierCoefficient[1] = vHat;
    fourierCoefficient[2] = wHat;


    return fourierCoefficient;
  }
}

/** \brief Generates a velocity field from Fourier-modes using FFTW
 * The disturbances are generated and filtered to a coarser resolution
 *
 * \author george
 * \date 24.01.2011
 *
 * (u,v,w)PhysFieldCoarse: real velocity in physical space
 *
 * \param uPhysFieldCoarse pointer to a MFloat for u
 * \param vPhysFieldCoarse pointer to a MFloat for v
 * \param wPhysFieldCoarse pointer to a MFloat for w
 *
 * ATTENTION: This version is quite old and should be thoroughly before use!
 */
inline void initFftFilter(MFloat* uPhysFieldCoarse,
                          MFloat* vPhysFieldCoarse,
                          MFloat* wPhysFieldCoarse,
                          MInt lx,
                          MInt ly,
                          MInt lz,
                          MInt lxCoarse,
                          MInt lyCoarse,
                          MInt lzCoarse,
                          MInt noPeakModes,
                          const MFloat Ma) {
  TRACE();


  MFloat waveVector[3], k0;

  std::complex<MFloat>* fourierCoefficient;

  fftw_complex *uHatField, *vHatField, *wHatField;
  fftw_complex *uPhysField, *vPhysField, *wPhysField;

  fftw_plan planU, planV, planW;

  MInt xRatio = lx / lxCoarse;
  MInt yRatio = ly / lyCoarse;
  MInt zRatio = lz / lzCoarse;

  if(lx % 2 != 0 || ly % 2 != 0 || lz % 2 != 0) {
    mTerm(1, AT_, " FFTInit: domainsize must NOT be an odd number! ");
  }

  m_log << " --- initializing FFTW --- " << std::endl;
  m_log << " domain size = " << lx << "x" << ly << "x" << lz << std::endl;

  // allocate field of velocities (is deleted at the end of the method)
  uPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  vPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  wPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

  // allocate field of Fourier coefficients (is deleted after the transformation)
  uHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  vHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  wHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

  // create plans for FFTW
  planU = fftw_plan_dft_3d(lx, ly, lz, uHatField, uPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planV = fftw_plan_dft_3d(lx, ly, lz, vHatField, vPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planW = fftw_plan_dft_3d(lx, ly, lz, wHatField, wPhysField, FFTW_BACKWARD, FFTW_MEASURE);

  // reset coefficients and velocities
  for(MInt p = 0; p < lx; p++) {
    for(MInt q = 0; q < ly; q++) {
      for(MInt r = 0; r < lz; r++) {
        uHatField[r + lz * (q + ly * p)][0] = 0.0;
        uHatField[r + lz * (q + ly * p)][1] = 0.0;
        uPhysField[r + lz * (q + ly * p)][0] = 0.0;
        uPhysField[r + lz * (q + ly * p)][1] = 0.0;

        vHatField[r + lz * (q + ly * p)][0] = 0.0;
        vHatField[r + lz * (q + ly * p)][1] = 0.0;
        vPhysField[r + lz * (q + ly * p)][0] = 0.0;
        vPhysField[r + lz * (q + ly * p)][1] = 0.0;

        wHatField[r + lz * (q + ly * p)][0] = 0.0;
        wHatField[r + lz * (q + ly * p)][1] = 0.0;
        wPhysField[r + lz * (q + ly * p)][0] = 0.0;
        wPhysField[r + lz * (q + ly * p)][1] = 0.0;
      }
    }
  }

  // Set coefficients:
  // FFTW stores the coefficients for positive wavenumbers in the first half of the array,
  // and those for negative wavenumbers in reverse order in the second half.
  // [0, 1, ... , N/2-1, N/2, ... , N-1]
  //  - the entry for zero-wavenumber is at position 0
  //  - the k-th entry and the (N-k)th entry correspond to wavenumbers with opposite sign
  //  - the entry at position N/2 corresponds to the Nyquist wavenumber and appears only once

  // peak wave number of energy spectrum
  k0 = 2.0 * PI / (lx / noPeakModes);

  for(MInt p = 0; p <= lx / 2; p++) {
    for(MInt q = 0; q <= ly / 2; q++) {
      for(MInt r = 0; r <= lz / 2; r++) {
        // wave-vector: k(p,q,r) = (2 \pi p / lx, 2 \pi q / ly, 2 \pi r / lz)
        waveVector[0] = (p)*2.0 * PI / lx;
        waveVector[1] = (q)*2.0 * PI / ly;
        waveVector[2] = (r)*2.0 * PI / lz;

        fourierCoefficient = getFourierCoefficient(waveVector, k0, Ma);

        // 1. Positive frequencies:
        uHatField[r + lz * (q + ly * p)][0] = std::real(fourierCoefficient[0]);
        uHatField[r + lz * (q + ly * p)][1] = std::imag(fourierCoefficient[0]);

        vHatField[r + lz * (q + ly * p)][0] = std::real(fourierCoefficient[1]);
        vHatField[r + lz * (q + ly * p)][1] = std::imag(fourierCoefficient[1]);

        wHatField[r + lz * (q + ly * p)][0] = std::real(fourierCoefficient[2]);
        wHatField[r + lz * (q + ly * p)][1] = std::imag(fourierCoefficient[2]);

        // 2. Negative frequencies:
        if(p > 1 && q > 1 && r > 1) {
          if(p < lx / 2 && q < ly / 2 && r < lz / 2) {
            // since the physical velocity field is real, the coefficients for negative frequencies
            // are the std::complex conjugate of those for positive frequencies
            uHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = uHatField[r + lz * (q + ly * p)][0];
            uHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -uHatField[r + lz * (q + ly * p)][1];

            vHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = vHatField[r + lz * (q + ly * p)][0];
            vHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -vHatField[r + lz * (q + ly * p)][1];

            wHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = wHatField[r + lz * (q + ly * p)][0];
            wHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -wHatField[r + lz * (q + ly * p)][1];
          }
        }
      }
    }
  }

  // Do Fourier transform (backward, see plan definition)
  // Definition in one dimension:
  // u(x) = \sum_{j=0}^{lx-1} \hat{u}_j exp(i 2 \pi j x / lx)
  fftw_execute(planU);
  fftw_execute(planV);
  fftw_execute(planW);

  // normalize (this preserves the norm of the basis functions)
  for(MInt p = 0; p < lx; p++) {
    for(MInt q = 0; q < ly; q++) {
      for(MInt r = 0; r < lz; r++) {
        uPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));
        vPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));
        wPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));

        uPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
        vPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
        wPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
      }
    }
  }

  // // test:
  // // cos(x) in a cube with 32x32x32 cells
  // uHatField[(0+lz*(0+ly*2)=2048][0]=90.5;
  // uHatField[(0+lz*(0+ly*30)=30720][0]=90.5;

  // reset coarse fields
  for(MInt p = 0; p < lxCoarse; p++) {
    for(MInt q = 0; q < lyCoarse; q++) {
      for(MInt r = 0; r < lzCoarse; r++) {
        uPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] = F0;
        vPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] = F0;
        wPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] = F0;
      }
    }
  }

  // transfer physField to physFieldCoarse via arithmetic averaging
  for(MInt p = 0; p < lxCoarse; p++) {
    for(MInt q = 0; q < lyCoarse; q++) {
      for(MInt r = 0; r < lzCoarse; r++) {
        for(MInt i = p * xRatio; i < p * zRatio + xRatio; i++) {
          for(MInt j = q * yRatio; j < q * zRatio + yRatio; j++) {
            for(MInt k = r * zRatio; k < r * zRatio + zRatio; k++) {
              uPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] += uPhysField[k + lz * (j + ly * i)][0];
              vPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] += vPhysField[k + lz * (j + ly * i)][0];
              wPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] += wPhysField[k + lz * (j + ly * i)][0];
            }
          }
        }
      }
    }
  }
  for(MInt p = 0; p < lxCoarse; p++) {
    for(MInt q = 0; q < lyCoarse; q++) {
      for(MInt r = 0; r < lzCoarse; r++) {
        uPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] /= (xRatio * yRatio * zRatio);
        vPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] /= (xRatio * yRatio * zRatio);
        wPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] /= (xRatio * yRatio * zRatio);
      }
    }
  }

  fftw_destroy_plan(planU);
  fftw_destroy_plan(planV);
  fftw_destroy_plan(planW);

  fftw_free(uHatField);
  fftw_free(vHatField);
  fftw_free(wHatField);

  fftw_free(uPhysField);
  fftw_free(vPhysField);
  fftw_free(wPhysField);
}

/** \brief Generates a velocity field from Fourier-modes using FFTW
 *
 * \author george
 * \date 09.12.2011
 *
 * (u,v,w) PhysField: complex velocity in physical space
 *
 * \param uPhysField pointer to a fftw_complex for u
 * \param vPhysField pointer to a fftw_complex for v
 * \param wPhysField pointer to a fftw_complex for w
 *
 * ATTENTION: This version is quite old and should be thoroughly before use!
 */
inline void initFft(fftw_complex* uPhysField,
                    fftw_complex* vPhysField,
                    fftw_complex* wPhysField,
                    MInt lx,
                    MInt ly,
                    MInt lz,
                    MInt noPeakModes,
                    const MFloat Ma) {
  TRACE();

  MFloat waveVector[3], k0;

  std::complex<MFloat>* fourierCoefficient;

  fftw_complex *uHatField, *vHatField, *wHatField;

  fftw_plan planU, planV, planW;

  if(lx % 2 != 0 || ly % 2 != 0 || lz % 2 != 0) {
    std::stringstream errorMessage;
    errorMessage << " FFTInit: domainsize must NOT be an odd number! " << std::endl;
    mTerm(1, AT_, errorMessage.str());
  }

  m_log << " --- initializing FFTW --- " << std::endl;
  m_log << " domain size = " << lx << "x" << ly << "x" << lz << std::endl;

  // allocate field of Fourier coefficients (is deleted after the transformation)
  uHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  vHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  wHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

  // create plans for FFTW
  planU = fftw_plan_dft_3d(lx, ly, lz, uHatField, uPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planV = fftw_plan_dft_3d(lx, ly, lz, vHatField, vPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planW = fftw_plan_dft_3d(lx, ly, lz, wHatField, wPhysField, FFTW_BACKWARD, FFTW_MEASURE);

  // reset coefficients and velocities
  for(MInt p = 0; p < lx; p++) {
    for(MInt q = 0; q < ly; q++) {
      for(MInt r = 0; r < lz; r++) {
        uHatField[r + lz * (q + ly * p)][0] = 0.0;
        uHatField[r + lz * (q + ly * p)][1] = 0.0;
        uPhysField[r + lz * (q + ly * p)][0] = 0.0;
        uPhysField[r + lz * (q + ly * p)][1] = 0.0;

        vHatField[r + lz * (q + ly * p)][0] = 0.0;
        vHatField[r + lz * (q + ly * p)][1] = 0.0;
        vPhysField[r + lz * (q + ly * p)][0] = 0.0;
        vPhysField[r + lz * (q + ly * p)][1] = 0.0;

        wHatField[r + lz * (q + ly * p)][0] = 0.0;
        wHatField[r + lz * (q + ly * p)][1] = 0.0;
        wPhysField[r + lz * (q + ly * p)][0] = 0.0;
        wPhysField[r + lz * (q + ly * p)][1] = 0.0;
      }
    }
  }

  // Set coefficients:
  // FFTW stores the coefficients for positive wavenumbers in the first half of the array,
  // and those for negative wavenumbers in reverse order in the second half.
  // [0, 1, ... , N/2-1, N/2, ... , N-1]
  //  - the entry for zero-wavenumber is at position 0
  //  - the k-th entry and the (N-k)th entry correspond to wavenumbers with opposite sign
  //  - the entry at position N/2 corresponds to the Nyquist wavenumber and appears only once

  // peak wave number of energy spectrum
  k0 = 2.0 * PI / (lx / noPeakModes);

  for(MInt p = 0; p <= lx / 2; p++) {
    for(MInt q = 0; q <= ly / 2; q++) {
      for(MInt r = 0; r <= lz / 2; r++) {
        // wave-vector: k(p,q,r) = (2 \pi p / lx, 2 \pi q / ly, 2 \pi r / lz)
        waveVector[0] = (p)*2.0 * PI / lx;
        waveVector[1] = (q)*2.0 * PI / ly;
        waveVector[2] = (r)*2.0 * PI / lz;

        fourierCoefficient = getFourierCoefficient(waveVector, k0, Ma);

        // 1. Positive frequencies:
        uHatField[r + lz * (q + ly * p)][0] = std::real(fourierCoefficient[0]);
        uHatField[r + lz * (q + ly * p)][1] = std::imag(fourierCoefficient[0]);

        vHatField[r + lz * (q + ly * p)][0] = std::real(fourierCoefficient[1]);
        vHatField[r + lz * (q + ly * p)][1] = std::imag(fourierCoefficient[1]);

        wHatField[r + lz * (q + ly * p)][0] = std::real(fourierCoefficient[2]);
        wHatField[r + lz * (q + ly * p)][1] = std::imag(fourierCoefficient[2]);

        // 2. Negative frequencies:
        if(p > 1 && q > 1 && r > 1) {
          if(p < lx / 2 && q < ly / 2 && r < lz / 2) {
            // since the physical velocity field is real, the coefficients for negative frequencies
            // are the std::complex conjugate of those for positive frequencies
            uHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = uHatField[r + lz * (q + ly * p)][0];
            uHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -uHatField[r + lz * (q + ly * p)][1];

            vHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = vHatField[r + lz * (q + ly * p)][0];
            vHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -vHatField[r + lz * (q + ly * p)][1];

            wHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = wHatField[r + lz * (q + ly * p)][0];
            wHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -wHatField[r + lz * (q + ly * p)][1];
          }
        }
      }
    }
  }

  // // test:
  // // cos(x) in a cube with 32x32x32 cells
  // uHatField[(0+lz*(0+ly*2)=2048][0]=90.5;
  // uHatField[(0+lz*(0+ly*30)=30720][0]=90.5;

  // Do Fourier transform (backward, see plan definition)
  // Definition in one dimension:
  // u(x) = \sum_{j=0}^{lx-1} \hat{u}_j exp(i 2 \pi j x / lx)
  fftw_execute(planU);
  fftw_execute(planV);
  fftw_execute(planW);

  // normalize (this preserves the norm of the basis functions)
  for(MInt p = 0; p < lx; p++) {
    for(MInt q = 0; q < ly; q++) {
      for(MInt r = 0; r < lz; r++) {
        uPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));
        vPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));
        wPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));

        uPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
        vPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
        wPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
      }
    }
  }

  fftw_destroy_plan(planU);
  fftw_destroy_plan(planV);
  fftw_destroy_plan(planW);
  fftw_free(uHatField);
  fftw_free(vHatField);
  fftw_free(wHatField);
}


/** \brief Compute energy spectrum on unity cube
 *
 * \author Lennart Schneiders
 *
 * discrete spectrum is computed on unity cube, no spatial scaling required
 * velocity field is computed for u_rms = 1, hence uPhysField subsequently has to be scaled by the magnitude of the
 * fluctuations, e.g. uPhysField *= m_UInfinity
 *
 * (u,v,w) PhysField: complex velocity in physical space
 *
 * \param uPhysField pointer to a fftw_complex for u
 * \param vPhysField pointer to a fftw_complex for v
 * \param wPhysField pointer to a fftw_complex for w
 */
inline void computeEnergySpectrum(MFloatScratchSpace& velocity, MIntScratchSpace& indices, const ptrdiff_t howmany,
                                  MInt locSize, MInt nx, MInt ny, MInt nz, MFloat viscosity, const MPI_Comm comm) {
  DEBUG("computeEnergySpectrum entry", MAIA_DEBUG_TRACE_IN);
  TRACE();

  const MInt size = nx * ny * nz;
  const MFloat fsize = (MFloat)size;
  const ptrdiff_t rank = 3;
  if(nx % 2 != 0 || ny % 2 != 0 || nz % 2 != 0) {
    std::stringstream errorMessage;
    errorMessage << " FFTInit: domainsize must NOT be an odd number! " << nx << " " << ny << " " << nz << std::endl;
    mTerm(1, AT_, errorMessage.str());
  }
  if(howmany < 3) mTerm(1, AT_, "expecting at least 3D velocity field");

  m_log << " --- initializing FFTW --- " << std::endl;
  m_log << " domain size = " << nx << "x" << ny << "x" << nz << " (" << howmany << ")" << std::endl;

  // caution: nz must be evenly divisble by number of ranks!
  MInt domainId;
  MInt noDomains;
  MPI_Comm_rank(comm, &domainId);
  MPI_Comm_size(comm, &noDomains);
  MPI_Comm MPI_COMM_FFTW;
  MInt maxRank = 1;
  maxRank = mMin(nz, noDomains);
  while(nz % maxRank != 0) {
    maxRank--;
  }
  MInt color = (domainId < maxRank) ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(comm, color, domainId, &MPI_COMM_FFTW, AT_, "MPI_COMM_FFTW");
  m_log << "no ranks for fft: " << maxRank << std::endl;

  MFloat urms0 = F0;
  for(MInt i = 0; i < locSize; i++) {
    urms0 += POW2(velocity(i, 0)) + POW2(velocity(i, 1)) + POW2(velocity(i, 2));
  }
  MPI_Allreduce(MPI_IN_PLACE, &urms0, 1, MPI_DOUBLE, MPI_SUM, comm, AT_, "MPI_IN_PLACE", "urms0");
  urms0 = sqrt(F1B3 * urms0 / fsize);

  MInt fftLocSize = 0;
  ptrdiff_t alloc_local = 0, local_n0 = 0, local_0_start = 0;
  fftw_complex* hatField = nullptr;
  fftw_plan plan;

  const MFloat time0 = MPI_Wtime();

  const ptrdiff_t n[3] = {nx, ny, nz};
#ifndef MAIA_WINDOWS
  alloc_local = (domainId < maxRank) ? fftw_mpi_local_size_many(rank, n, howmany, FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_FFTW,
                                                                &local_n0, &local_0_start)
                                     : 1;
#endif

  ScratchSpace<fftw_complex> hatFieldMem(alloc_local, AT_, "hatFieldMem");
  if(domainId < maxRank) {
    hatField = &hatFieldMem[0];
  } else {
    alloc_local = 0;
    local_n0 = 0;
    local_0_start = nx;
  }

  ScratchSpace<MLong> offsets(noDomains + 1, AT_, "offsets");
  ScratchSpace<MInt> noRecv(noDomains, AT_, "noRecv");
  ScratchSpace<MInt> noSend(noDomains, AT_, "noSend");
  noRecv.fill(0);
  noSend.fill(0);
  MLong locOffset = (domainId < maxRank) ? ((MInt)local_0_start) * ny * nz : size;

  MPI_Allgather(&locOffset, 1, MPI_LONG, &offsets[0], 1, MPI_LONG, comm, AT_, "locOffset", "offsets[0]");
  offsets(noDomains) = size;
  if(domainId < maxRank) {
    if(offsets(domainId) != locOffset || offsets(domainId + 1) != ((MInt)(local_0_start + local_n0)) * ny * nz)
      mTerm(1, AT_, "wrong size 0");
  }
  for(MInt i = 0; i < locSize; i++) {
    MInt j = indices(i);
    MInt nghbrDomain = mMin(maxRank - 1, j / (size / maxRank));
    while(j < offsets(nghbrDomain) || j >= offsets(nghbrDomain + 1)) {
      if(j < offsets(nghbrDomain)) nghbrDomain--;
      if(j >= offsets(nghbrDomain + 1)) nghbrDomain++;
    }
    if(nghbrDomain >= maxRank) mTerm(1, AT_, "wrong domain");
    noSend(nghbrDomain)++;
  }
  MPI_Alltoall(&noSend[0], 1, MPI_INT, &noRecv[0], 1, MPI_INT, comm, AT_, "noSend[0]", "noRecv[0]");
  MInt sendSize = 0;
  MInt recvSize = 0;
  ScratchSpace<MInt> recvOffsets(noDomains + 1, AT_, "recvOffsets");
  ScratchSpace<MInt> sendOffsets(noDomains + 1, AT_, "sendOffsets");
  for(MInt i = 0; i < noDomains; i++) {
    recvOffsets(i) = recvSize;
    sendOffsets(i) = sendSize;
    sendSize += noSend(i);
    recvSize += noRecv(i);
  }
  recvOffsets(noDomains) = recvSize;
  sendOffsets(noDomains) = sendSize;
  ScratchSpace<MFloat> recvData(mMax(1, recvSize), howmany, AT_, "recvData");
  ScratchSpace<MFloat> sendData(mMax(1, sendSize), howmany, AT_, "sendData");
  ScratchSpace<MInt> recvIndices(mMax(1, recvSize), AT_, "recvIndices");
  ScratchSpace<MInt> sendIndices(mMax(1, sendSize), AT_, "sendIndices");
  noSend.fill(0);
  for(MInt i = 0; i < locSize; i++) {
    MInt j = indices(i);
    MInt nghbrDomain = mMin(maxRank - 1, j / (size / maxRank));
    while(j < offsets(nghbrDomain) || j >= offsets(nghbrDomain + 1)) {
      if(j < offsets(nghbrDomain)) nghbrDomain--;
      if(j >= offsets(nghbrDomain + 1)) nghbrDomain++;
    }
    if(nghbrDomain >= maxRank) mTerm(1, AT_, "wrong domain");
    if(sendOffsets(nghbrDomain) + noSend(nghbrDomain) >= sendOffsets(nghbrDomain + 1)) mTerm(1, AT_, "wrong size");
    for(MInt k = 0; k < howmany; k++) {
      sendData(sendOffsets(nghbrDomain) + noSend(nghbrDomain), k) = velocity(i, k);
    }
    sendIndices(sendOffsets(nghbrDomain) + noSend(nghbrDomain)) = j;
    if(j < offsets(nghbrDomain) || j >= offsets(nghbrDomain + 1)) mTerm(1, AT_, "wrong size 01");
    noSend(nghbrDomain)++;
  }
  for(MInt i = 0; i < noDomains; i++) {
    if(noSend(i) != sendOffsets(i + 1) - sendOffsets(i)) mTerm(1, AT_, "wrong size 1");
    if(noRecv(i) != recvOffsets(i + 1) - recvOffsets(i)) mTerm(1, AT_, "wrong size 1");
  }

  ScratchSpace<MPI_Request> sendReq(noDomains, AT_, "sendReq");
  sendReq.fill(MPI_REQUEST_NULL);
  MInt cnt = 0;
  for(MInt i = 0; i < noDomains; i++) {
    if(noSend(i) == 0) continue;
    MPI_Issend(&sendData[howmany * sendOffsets(i)], howmany * noSend(i), MPI_DOUBLE, i, 123, comm, &sendReq[cnt], AT_,
               "sendData[howmany * sendOffsets(i)]");
    cnt++;
  }
  for(MInt i = 0; i < noDomains; i++) {
    if(noRecv(i) == 0) continue;
    MPI_Recv(&recvData[howmany * recvOffsets(i)], howmany * noRecv(i), MPI_DOUBLE, i, 123, comm, MPI_STATUS_IGNORE, AT_,
             "recvData[howmany * recvOffsets(i)]");
  }
  if(cnt > 0) MPI_Waitall(cnt, &sendReq[0], MPI_STATUSES_IGNORE, AT_);
  sendReq.fill(MPI_REQUEST_NULL);
  cnt = 0;
  for(MInt i = 0; i < noDomains; i++) {
    if(noSend(i) == 0) continue;
    MPI_Issend(&sendIndices[sendOffsets(i)], noSend(i), MPI_INT, i, 124, comm, &sendReq[cnt], AT_,
               "sendIndices[sendOffsets(i)]");
    cnt++;
  }
  for(MInt i = 0; i < noDomains; i++) {
    if(noRecv(i) == 0) continue;
    MPI_Recv(&recvIndices[recvOffsets(i)], noRecv(i), MPI_INT, i, 124, comm, MPI_STATUS_IGNORE, AT_,
             "recvIndices[recvOffsets(i)]");
  }
  if(cnt > 0) MPI_Waitall(cnt, &sendReq[0], MPI_STATUSES_IGNORE, AT_);

  if(domainId < maxRank) {
    for(MInt i = 0; i < recvSize; i++) {
      MInt j = recvIndices(i);
      if(j < ((MInt)local_0_start) * ny * nz || j >= ((MInt)(local_0_start + local_n0)) * ny * nz)
        mTerm(1, AT_, "wrong size 2");
      MInt pos = j - (((MInt)local_0_start) * ny * nz);
      for(MInt k = 0; k < howmany; k++) {
        hatField[howmany * pos + k][0] = recvData(i, k); // Real part of fourier transformed velocity
        hatField[howmany * pos + k][1] = F0;             // Imaginary part of fourier transformed velocity
      }
    }
    fftLocSize = recvSize;

    MFloat couplcheck = F0;
    for(MInt i = 0; i < fftLocSize; i++) {
      for(MInt k = 0; k < 3; k++) {
        couplcheck -= hatField[howmany * i + k][0] * hatField[howmany * i + 3 + k][0];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &couplcheck, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE", "couplcheck");
    if(domainId == 0)
      std::cerr << "couplecheck fft " << couplcheck << " " << couplcheck * POW3(0.25 / 64.0) << std::endl;
  }

  // Perform Fourier transform
  if(domainId < maxRank) {
#ifndef MAIA_WINDOWS
    plan = fftw_mpi_plan_many_dft(rank, n, howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, hatField, hatField,
                                  MPI_COMM_FFTW, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
#endif
    for(MInt i = 0; i < fftLocSize; i++) {
      for(MInt j = 0; j < howmany; j++) {
        hatField[howmany * i + j][0] /= fsize;
        hatField[howmany * i + j][1] /= fsize;
      }
    }
  }

  const MFloat time1 = MPI_Wtime();

  m_log << "parallel fftw time " << time1 - time0 << std::endl;

  // Set coefficients:
  // FFTW stores the coefficients for positive wavenumbers in the first half of the array,
  // and those for negative wavenumbers in reverse order in the second half.
  // [0, 1, ... , N/2-1, N/2, ... , N-1]
  //  - the entry for zero-wavenumber is at position 0
  //  - the k-th entry and the (N-k)th entry correspond to wavenumbers with opposite sign
  //  - the entry at position N/2 corresponds to the Nyquist wavenumber and appears only once

  static const MInt referenceCubeSize =
      ((Context::propertyExists("referenceCubeSize")) ? Context::getBasicProperty<MInt>("referenceCubeSize", AT_) : 1);

  static const MInt kMinSpec = ((Context::propertyExists("kMin")) ? Context::getBasicProperty<MInt>("kMin", AT_) : 0);

  static const MInt kMaxSpec = ((Context::propertyExists("kMax")) ? Context::getBasicProperty<MInt>("kMax", AT_) : 0);

  const MFloat k0 = F2 * PI / referenceCubeSize;
  MInt noBins = 0;
  MFloat kMin = 0.0;
  MFloat kMax = 0.0;
  MFloat deltaK = 0.0;
  if(kMinSpec > 0 && kMaxSpec > 0) {
    noBins = kMaxSpec;
    kMin = static_cast<float>(kMinSpec) * k0;
    kMax = static_cast<float>(kMaxSpec) * k0;
    deltaK = (kMax - kMin) / ((MFloat)(kMaxSpec - 1));
  } else {
    noBins = mMax(nx, mMax(ny, nz)) / 2;
    kMin = F1 * k0;
    kMax = noBins * k0;
    deltaK = (kMax - kMin) / ((MFloat)(noBins - 1));
  }

  // const MInt noBins = 200;
  // const MFloat kMax = 200.0*k0;
  ScratchSpace<MFloat> spectrumBnd(noBins + 1, AT_, "spectrumBnd");
  ScratchSpace<MFloat> spectrum(noBins, AT_, "spectrum");
  ScratchSpace<MFloat> coupling(noBins, AT_, "coupling");
  ScratchSpace<MFloat> transfer(noBins, AT_, "transfer");
  ScratchSpace<MFloat> transfer2(noBins, AT_, "transfer2");
  ScratchSpace<MFloat> rate(noBins, AT_, "rate");
  ScratchSpace<MFloat> flux(noBins, 2, AT_, "flux");
  ScratchSpace<MFloat> flux2(noBins, 2, AT_, "flux2");
  ScratchSpace<MFloat> spectrumSum(noBins, AT_, "spectrumSum");
  ScratchSpace<MFloat> spectrumSum2(noBins, AT_, "spectrumSum2");
  spectrumBnd[0] = kMin - F1B2 * deltaK;
  const MFloat spectrumBnd00 = kMin - F3B4 * deltaK;
  for(MInt i = 0; i < noBins; i++) {
    spectrumBnd[i + 1] = spectrumBnd[i] + deltaK;
  }
  spectrum.fill(F0);
  coupling.fill(F0);
  transfer.fill(F0);
  transfer2.fill(F0);
  rate.fill(F0);
  flux.fill(F0);
  flux2.fill(F0);
  spectrumSum.fill(F0);
  spectrumSum2.fill(F0);

  // beware: very expensive routine
  static const MBool computeTransferRate =
      ((Context::propertyExists("computeTransferRate")) ? Context::getBasicProperty<MBool>("computeTransferRate", AT_)
                                                        : false);
  if(computeTransferRate) {
    MIntScratchSpace noDatPerDomain(noDomains, AT_, "noDatPerDomain");
    MIntScratchSpace dataOffsets(noDomains + 1, AT_, "dataOffsets");
    MPI_Allgather(&fftLocSize, 1, MPI_INT, &noDatPerDomain[0], 1, MPI_INT, comm, AT_, "fftLocSize",
                  "noDatPerDomain[0]");
    dataOffsets.fill(0);
    for(MInt d = 0; d < noDomains; d++) {
      dataOffsets(d + 1) = dataOffsets(d) + noDatPerDomain(d);
    }

    MPI_Request sendReqLoc = MPI_REQUEST_NULL;
    fftw_complex* uHatGlobal = fftw_alloc_complex(3 * size);
    fftw_complex* uHatGlobalTmp = nullptr;
    for(MInt i = 0; i < fftLocSize; i++) {
      for(MInt j = 0; j < 3; j++) {
        uHatGlobal[3 * i + j][0] = hatField[howmany * i + j][0];
        uHatGlobal[3 * i + j][1] = hatField[howmany * i + j][1];
      }
    }
    if(fftLocSize > 0) {
      MPI_Issend(&uHatGlobal[0], 6 * fftLocSize, MPI_DOUBLE, 0, 18, comm, &sendReqLoc, AT_, "uHatGlobal[0]");
    }
    if(domainId == 0) {
      uHatGlobalTmp = fftw_alloc_complex(3 * size);
      for(MInt d = 0; d < noDomains; d++) {
        if(noDatPerDomain[d] > 0) {
          MPI_Recv(&(uHatGlobalTmp[3 * dataOffsets[d]][0]), 6 * noDatPerDomain[d], MPI_DOUBLE, d, 18, comm,
                   MPI_STATUSES_IGNORE, AT_, "(uHatGlobalTmp[3 * dataOffsets[d]][0])");
        }
      }
    }
    if(fftLocSize > 0) {
      MPI_Wait(&sendReqLoc, MPI_STATUSES_IGNORE, AT_);
    }

    if(domainId == 0) {
      // re-sort by ascending wave numbers
      for(MInt i0 = 0; i0 < nx; i0++) {
        for(MInt i1 = 0; i1 < ny; i1++) {
          for(MInt i2 = 0; i2 < nz; i2++) {
            MInt j0 = (i0 > nx / 2) ? i0 - nx : i0;
            MInt j1 = (i1 > ny / 2) ? i1 - ny : i1;
            MInt j2 = (i2 > nz / 2) ? i2 - nz : i2;
            j0 += nx / 2 - 1;
            j1 += nx / 2 - 1;
            j2 += nx / 2 - 1;
            MInt pos0 = getGlobalPosFFTW(i0, i1, i2, ny, nz);
            MInt pos1 = getGlobalPosFFTW(j0, j1, j2, ny, nz);
            for(MInt j = 0; j < 3; j++) {
              uHatGlobal[3 * pos1 + j][0] = uHatGlobalTmp[3 * pos0 + j][0];
              uHatGlobal[3 * pos1 + j][1] = uHatGlobalTmp[3 * pos0 + j][1];
            }
          }
        }
      }
      fftw_free(uHatGlobalTmp);
    }

    MPI_Bcast(&uHatGlobal[0], 6 * size, MPI_DOUBLE, 0, comm, AT_, "uHatGlobal");


    // MInt domainSize = size/noDomains;
    // MInt posMin = domainSize*domainId;
    // MInt posMax = mMin(size,domainSize*(domainId+1));
    // cerr << domainId << ": pos " << posMin << " " << posMax << std::endl;

    MInt noWaveNumbers = 0;
    MLong weightSum = 0;
    for(MInt pos = 0; pos < size; pos++) {
      MInt i2 = pos % nz;
      MInt i1 = (pos / nz) % ny;
      MInt i0 = pos / (ny * nz);
      if(pos != i2 + nz * (i1 + ny * (i0))) mTerm(1, AT_, "index mismatch");
      MFloat k[3];
      MInt j0 = i0 - nx / 2 + 1;
      MInt j1 = i1 - ny / 2 + 1;
      MInt j2 = i2 - nz / 2 + 1;
      k[0] = ((MFloat)j0) * k0;
      k[1] = ((MFloat)j1) * k0;
      k[2] = ((MFloat)j2) * k0;

      MFloat kAbs = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
      if(kAbs > spectrumBnd[noBins] + 1e-12) continue;
      if(kAbs < spectrumBnd[0] - 1e-12) continue;

      MInt p0Min = std::max(0, i0 - nx / 2);
      MInt p0Max = std::min(nx, i0 + nx / 2);
      MInt p1Min = std::max(0, i1 - ny / 2);
      MInt p1Max = std::min(ny, i1 + ny / 2);
      MInt p2Min = std::max(0, i2 - nz / 2);
      MInt p2Max = std::min(nz, i2 + nz / 2);
      MLong weight = ((p0Max - p0Min) * (p1Max - p1Min) * (p2Max - p2Min));

      weightSum += weight;
      noWaveNumbers++;
    }
    // MInt domainSize = noWaveNumbers/noDomains;
    MLong noDomainsLong = (MLong)noDomains;
    MLong domainIdLong = (MLong)domainId;
    MLong domainSize = weightSum / noDomainsLong;
    MInt posMin = 0;
    MInt posMax = size;
    noWaveNumbers = 0;
    weightSum = 0;
    for(MInt pos = 0; pos < size; pos++) {
      MInt i2 = pos % nz;
      MInt i1 = (pos / nz) % ny;
      MInt i0 = pos / (ny * nz);
      if(pos != i2 + nz * (i1 + ny * (i0))) mTerm(1, AT_, "index mismatch");
      MFloat k[3];
      MInt j0 = i0 - nx / 2 + 1;
      MInt j1 = i1 - ny / 2 + 1;
      MInt j2 = i2 - nz / 2 + 1;
      k[0] = ((MFloat)j0) * k0;
      k[1] = ((MFloat)j1) * k0;
      k[2] = ((MFloat)j2) * k0;

      MFloat kAbs = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
      if(kAbs > spectrumBnd[noBins] + 1e-12) continue;
      if(kAbs < spectrumBnd[0] - 1e-12) continue;

      MInt p0Min = std::max(0, i0 - nx / 2);
      MInt p0Max = std::min(nx, i0 + nx / 2);
      MInt p1Min = std::max(0, i1 - ny / 2);
      MInt p1Max = std::min(ny, i1 + ny / 2);
      MInt p2Min = std::max(0, i2 - nz / 2);
      MInt p2Max = std::min(nz, i2 + nz / 2);
      MLong weight = ((p0Max - p0Min) * (p1Max - p1Min) * (p2Max - p2Min));

      if(weightSum <= (domainSize * domainIdLong) && (weightSum + weight) > (domainSize * domainIdLong)) posMin = pos;
      if(weightSum <= (domainSize * (domainIdLong + 1)) && (weightSum + weight) > (domainSize * (domainIdLong + 1)))
        posMax = pos;

      weightSum += weight;
      noWaveNumbers++;
    }

    weightSum = 0;
    for(MInt pos = posMin; pos < posMax; pos++) {
      const MInt i2 = pos % nz;
      const MInt i1 = (pos / nz) % ny;
      const MInt i0 = pos / (ny * nz);

      const MInt j0 = i0 - nx / 2 + 1;
      const MInt j1 = i1 - ny / 2 + 1;
      const MInt j2 = i2 - nz / 2 + 1;
      const MFloat k[3] = {((MFloat)j0) * k0, ((MFloat)j1) * k0, ((MFloat)j2) * k0};
      const MFloat kAbs2 = mMax(1e-12, k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
      const MFloat kAbs = mMax(1e-12, sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]));

      if(kAbs > spectrumBnd[noBins] + 1e-12) continue;
      const MInt bin = ((MInt)((kAbs - spectrumBnd[0] + 1e-12) / deltaK));
      const MInt bin2 = ((MInt)((kAbs - spectrumBnd00 + 1e-12) / deltaK));

      const MFloat kjl[3][3] = {
          {1.0 - (k[0] * k[0] / kAbs2), 0.0 - (k[0] * k[1] / kAbs2), 0.0 - (k[0] * k[2] / kAbs2)},
          {0.0 - (k[1] * k[0] / kAbs2), 1.0 - (k[1] * k[1] / kAbs2), 0.0 - (k[1] * k[2] / kAbs2)},
          {0.0 - (k[2] * k[0] / kAbs2), 0.0 - (k[2] * k[1] / kAbs2), 1.0 - (k[2] * k[2] / kAbs2)}};

      const MInt p0Min = std::max(0, i0 - nx / 2);
      const MInt p0Max = std::min(nx, i0 + nx / 2);
      const MInt p1Min = std::max(0, i1 - ny / 2);
      const MInt p1Max = std::min(ny, i1 + ny / 2);
      const MInt p2Min = std::max(0, i2 - nz / 2);
      const MInt p2Max = std::min(nz, i2 + nz / 2);

      MFloat trans = F0;

      for(MInt p0 = p0Min; p0 < p0Max; p0++) {
        for(MInt p1 = p1Min; p1 < p1Max; p1++) {
          for(MInt p2 = p2Min; p2 < p2Max; p2++) {
            const MInt pos1 = p2 + nz * (p1 + ny * (p0));
            const MInt s0 = i0 - p0 + nx / 2 - 1;
            const MInt s1 = i1 - p1 + ny / 2 - 1;
            const MInt s2 = i2 - p2 + nz / 2 - 1;
            const MInt pos2 = s2 + nz * (s1 + ny * (s0));

            for(MInt i = 0; i < 3; i++) {
              for(MInt j = 0; j < 3; j++) {
                for(MInt l = 0; l < 3; l++) {
                  trans += k[i] * kjl[j][l]
                           * triadImag(uHatGlobal[3 * pos1 + l], uHatGlobal[3 * pos2 + i], uHatGlobal[3 * pos + j]);
                }
              }
            }
          }
        }
      }

      if(bin > -1 && bin < noBins) transfer(bin) += trans;
      if(bin2 > -1 && bin2 < noBins) transfer2(bin2) += trans;
    }

    MPI_Allreduce(MPI_IN_PLACE, &transfer[0], noBins, MPI_DOUBLE, MPI_SUM, comm, AT_, "MPI_IN_PLACE", "transfer[0]");
    MPI_Allreduce(MPI_IN_PLACE, &transfer2[0], noBins, MPI_DOUBLE, MPI_SUM, comm, AT_, "MPI_IN_PLACE", "transfer2[0]");
    fftw_free(uHatGlobal);
  }

  const MFloat time2 = MPI_Wtime();

  if(domainId < maxRank) {
    for(MInt i0 = (MInt)local_0_start; i0 < (MInt)(local_0_start + local_n0); i0++) {
      for(MInt i1 = 0; i1 < ny; i1++) {
        for(MInt i2 = 0; i2 < nz; i2++) {
          MFloat k[3];
          MInt j0 = (i0 > nx / 2) ? i0 - nx : i0;
          MInt j1 = (i1 > ny / 2) ? i1 - ny : i1;
          MInt j2 = (i2 > nz / 2) ? i2 - nz : i2;
          k[0] = ((MFloat)j0) * k0;
          k[1] = ((MFloat)j1) * k0;
          k[2] = ((MFloat)j2) * k0;

          MInt pos = i2 + nz * (i1 + ny * (i0 - ((MInt)local_0_start)));
          if(pos < 0 || pos >= fftLocSize) mTerm(1, AT_, "wrong size 3");

          MFloat kAbs = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
          MFloat eng = F0;
          MFloat coupl = F0;
          MFloat dedt = F0;
          for(MInt q = 0; q < 3; q++) {
            eng += F1B2 * (POW2(hatField[howmany * pos + q][0]) + POW2(hatField[howmany * pos + q][1]));
          }

          // if not only velocity field is passed as argument, also coupling rate Psi(k) and change of Energy dE/dt is
          // computed
          if(howmany > 3) {
            for(MInt q = 0; q < 3; q++) {
              coupl += (hatField[howmany * pos + q][0] * hatField[howmany * pos + 3 + q][0]
                        + hatField[howmany * pos + q][1] * hatField[howmany * pos + 3 + q][1]);
              dedt += (hatField[howmany * pos + 6 + q][0] * hatField[howmany * pos + q][0]
                       + hatField[howmany * pos + 6 + q][1] * hatField[howmany * pos + q][1]);
            }
          }
          if(kAbs > spectrumBnd[noBins] + 1e-12) continue;

          MInt bin = ((MInt)((kAbs - spectrumBnd[0] + 1e-12) / deltaK));
          MInt bin2 = ((MInt)((kAbs - spectrumBnd00 + 1e-12) / deltaK));
          if(bin > -1 && bin < noBins) {
            spectrum(bin) += eng;
            coupling(bin) += coupl;
            rate(bin) += dedt;
            spectrumSum(bin) += F1;
          }
          if(bin2 > -1 && bin2 < noBins) {
            spectrumSum2(bin2) += F1;
          }
        }
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &spectrum[0], noBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE",
                  "spectrum[0]");
    MPI_Allreduce(MPI_IN_PLACE, &coupling[0], noBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE",
                  "coupling[0]");
    MPI_Allreduce(MPI_IN_PLACE, &rate[0], noBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE", "rate[0]");
    MPI_Allreduce(MPI_IN_PLACE, &spectrumSum[0], noBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE",
                  "spectrumSum[0]");
    MPI_Allreduce(MPI_IN_PLACE, &spectrumSum2[0], noBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, AT_, "MPI_IN_PLACE",
                  "spectrumSum2[0]");

    for(MInt i = 0; i < noBins; i++) {
      MFloat db = spectrumBnd00 - spectrumBnd[0];
      MFloat vol =
          F4B3 * PI * (POW3(spectrumBnd[i + 1]) - POW3(spectrumBnd[i])) / (spectrumBnd[i + 1] - spectrumBnd[i]);
      MFloat vol2 = F4B3 * PI * (POW3(spectrumBnd[i + 1] + db) - POW3(spectrumBnd[i] + db))
                    / (spectrumBnd[i + 1] - spectrumBnd[i]);
      spectrum(i) *= vol / (spectrumSum(i) * POW3(k0));
      coupling(i) *= vol / (spectrumSum(i) * POW3(k0));
      transfer(i) *= vol / (spectrumSum(i) * POW3(k0));
      transfer2(i) *= vol2 / (spectrumSum2(i) * POW3(k0));
      rate(i) *= vol / (spectrumSum(i) * POW3(k0));
    }

    for(MInt i = 0; i < noBins; i++) {
      flux(i, 0) = F0;
      flux(i, 1) = F0;
      for(MInt j = 0; j < i; j++) {
        flux(i, 0) += transfer(j) * (spectrumBnd[j + 1] - spectrumBnd[j]) / k0;
      }
      for(MInt j = i; j < noBins; j++) {
        flux(i, 1) += transfer(j) * (spectrumBnd[j + 1] - spectrumBnd[j]) / k0;
      }
      flux2(i, 0) = F0;
      flux2(i, 1) = F0;
      for(MInt j = 0; j < i; j++) {
        flux2(i, 0) += transfer2(j) * (spectrumBnd[j + 1] - spectrumBnd[j]) / k0;
      }
      for(MInt j = i; j < noBins; j++) {
        flux2(i, 1) += transfer2(j) * (spectrumBnd[j + 1] - spectrumBnd[j]) / k0;
      }
    }

    if(domainId == 0) {
      std::ofstream ofl;
      std::ofstream ofl2;
      ofl.open("./out/energySpectrum_00" + std::to_string(globalTimeStep), std::ios_base::out | std::ios_base::trunc);
      ofl2.open("./out/energySpectrum_coarse_00" + std::to_string(globalTimeStep),
                std::ios_base::out | std::ios_base::trunc);
      if(ofl.is_open() && ofl.good() && ofl2.is_open() && ofl2.good()) {
        MFloat energy = F0;
        MFloat dissip = F0;
        MFloat length = F0;
        MFloat dedt = F0;
        MFloat force = F0;
        MFloat urms = F0;
        MFloat transf = F0;
        MFloat transf2 = F0;
        for(MInt i = 0; i < noBins; i++) {
          urms += spectrum(i) * (spectrumBnd[i + 1] - spectrumBnd[i]);
        }
        urms = sqrt(F2B3 * urms); // warum hier 2/3 und nicht 1/3
        for(MInt i = 0; i < noBins; i++) {
          MFloat k = F1B2 * (spectrumBnd[i + 1] + spectrumBnd[i]);
          MFloat Ek = spectrum(i);
          if(std::isnan(Ek) || std::isnan(coupling(i)) || std::isnan(rate(i))) continue;
          energy += (spectrumBnd[i + 1] - spectrumBnd[i]) * Ek;
          dissip += F2 * viscosity * POW2(k) * (spectrumBnd[i + 1] - spectrumBnd[i]) * Ek;
          length += (F1B2 * PI / POW2(urms)) * (spectrumBnd[i + 1] - spectrumBnd[i]) * Ek / (k);
          force += (spectrumBnd[i + 1] - spectrumBnd[i]) * coupling(i);
          dedt += (spectrumBnd[i + 1] - spectrumBnd[i]) * rate(i);
          transf += (spectrumBnd[i + 1] - spectrumBnd[i]) * transfer(i);
          transf2 += (spectrumBnd[i + 1] - spectrumBnd[i]) * transfer2(i);
        }
        ofl << "# Energy spectrum at time step " << globalTimeStep << ", total energy Ek=" << energy
            << ", dissipation rate=" << dissip << ", u_rms=" << urms << " (" << urms0
            << "), integral length scale=" << length << ", Kolmogorov time scale=" << sqrt(viscosity / dissip)
            << ", coupling=" << force << ", dedt=" << dedt << ", transfer=" << transf << ", transfer2=" << transf2
            << std::endl;
        ofl << "# 1: wave number k" << std::endl;
        ofl << "# 2: energy E(k)" << std::endl;
        ofl << "# 3: viscous dissipation rate D(k)" << std::endl;
        ofl << "# 4: transfer rate T(k)" << std::endl;
        ofl << "# 5: coupling rate Psi(k)" << std::endl;
        ofl << "# 6: flux -F(k)" << std::endl;
        ofl << "# 7: flux F(k)" << std::endl;
        ofl << "# 8: dE/dt(k)" << std::endl;
        ofl << "# 9: alternative transfer rate T(k)" << std::endl;
        ofl << "# 10: alternative flux -F(k)" << std::endl;
        ofl << "# 11: alternative flux F(k)" << std::endl;
        ofl << "# 12: number of samples in corresponding bin" << std::endl;
        ofl << "# 13: lower boundary of corresponding bin" << std::endl;
        ofl << "# 14: upper boundary of corresponding bin" << std::endl;
        ofl << "# 15: number of samples other spectrum: " << std::endl;
        ofl << "# 16: 4 * PI * k^2" << std::endl;
        ofl << "# 17: 4/3 * PI * (upperBndry^3 - lowerBndry^3)" << std::endl;
        for(MInt i = 0; i < noBins; i++) {
          MFloat k = F1B2 * (spectrumBnd[i + 1] + spectrumBnd[i]);
          MFloat Ek = spectrum(i);
          MFloat coupl = coupling(i);
          MFloat trans = transfer(i);
          MFloat trans2 = transfer2(i);
          ofl << k / k0 << " " << Ek << " " << F2 * POW2(k) * viscosity * Ek << " " << trans << " " << coupl << " "
              << flux(i, 0) << " " << flux(i, 1) << " " << rate(i) << " " << trans2 << " " << flux2(i, 0) << " "
              << flux2(i, 1) << " " << spectrumSum(i) << " " << spectrumBnd[i] << " " << spectrumBnd[i + 1] << " "
              << spectrumSum2(i) << " " << F4 * PI * POW2(k) << " "
              << F4B3 * PI * (POW3(spectrumBnd[i + 1]) - POW3(spectrumBnd[i])) << std::endl;
        }
        for(MInt i = 0; i < noBins / 2; i++) {
          const MInt i0 = 2 * i;
          const MInt i1 = 2 * i + 1;
          const MFloat k = F1B2 * (spectrumBnd[i1 + 1] + spectrumBnd[i0]);
          MFloat vol0 =
              F4B3 * PI * (POW3(spectrumBnd[i0 + 1]) - POW3(spectrumBnd[i0])) / (spectrumBnd[i0 + 1] - spectrumBnd[i0]);
          MFloat vol1 =
              F4B3 * PI * (POW3(spectrumBnd[i1 + 1]) - POW3(spectrumBnd[i1])) / (spectrumBnd[i1 + 1] - spectrumBnd[i1]);
          MFloat vol =
              F4B3 * PI * (POW3(spectrumBnd[i1 + 1]) - POW3(spectrumBnd[i0])) / (spectrumBnd[i1 + 1] - spectrumBnd[i0]);
          MFloat f0 = spectrumSum(i0) * vol / (vol0 * (spectrumSum(i0) + spectrumSum(i1)));
          MFloat f1 = spectrumSum(i1) * vol / (vol1 * (spectrumSum(i0) + spectrumSum(i1)));

          ofl2 << k / k0 << " " << f1 * spectrum(i1) + f0 * spectrum(i0) << " "
               << F2 * POW2(k) * viscosity * (f1 * spectrum(i1) + f0 * spectrum(i0)) << " "
               << f1 * transfer(i1) + f0 * transfer(i0) << " " << f1 * coupling(i1) + f0 * coupling(i0) << " "
               << f1 * flux(i1, 0) + f0 * flux(i0, 0) << " " << f1 * flux(i1, 1) + f0 * flux(i0, 1) << " "
               << f1 * rate(i1) + f0 * rate(i0) << " " << f1 * transfer2(i1) + f0 * transfer2(i0) << " "
               << f1 * flux2(i1, 0) + f0 * flux2(i0, 0) << " " << f1 * flux2(i1, 1) + f0 * flux2(i0, 1) << " "
               << spectrumSum(i0) + spectrumSum(i1) << " " << spectrumBnd[i0] << " " << spectrumBnd[i1 + 1] << " "
               << spectrumSum2(i0) + spectrumSum2(i1) << " " << F4 * PI * POW2(k) << " "
               << F4B3 * PI * (POW3(spectrumBnd[i1 + 1]) - POW3(spectrumBnd[i0])) << std::endl;
        }
        ofl.close();
        ofl2.close();
      }
    }

    // cleanup
    // fftw_cleanup(); problems occur when calling clean up function!!!
  }

  const MFloat time3 = MPI_Wtime();

  if(domainId == 0)
    std::cerr << "fft time " << time1 - time0 << " " << time2 - time1 << " " << time3 - time2 << std::endl;

  m_log << "fft finished" << std::endl;

  DEBUG("computeEnergySpectrum return", MAIA_DEBUG_TRACE_OUT);
}

} // namespace math
} // namespace maia
#endif // MAIA_FFTW_H
