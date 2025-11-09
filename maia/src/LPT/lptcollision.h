// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_PARTICLECOLLISION_H
#define MAIA_PARTICLECOLLISION_H

#include <vector>
#include "GRID/cartesiangridproxy.h"
#include "INCLUDE/maiatypes.h"
#include "lptcollisiondata.h"
#include "lptlib.h"
template <MInt nDim>
class ParticleCollision {
 private:
  using GridProxy = typename maia::grid::Proxy<nDim>;
  friend class LPT<nDim>;
  using lptParticleI = typename maia::lpt::partListIteratorConst<nDim>;
  using lptEllipticleI = typename maia::lpt::ellipsListIteratorConst<nDim>;

 public:
  ParticleCollision(LPT<nDim>* lptSolver, MInt collisionModel, const MInt solverId, const MInt domainId,
                    const MInt noDomains, const MPI_Comm comm);

  void init();

  void detectPartColl(std::vector<LPTSpherical<nDim>>& partList, std::vector<LPTEllipsoidal<nDim>>& partListEllipsoid,
                      const MFloat timeStep, const MFloat time);

  void geometryInteraction();

  void writeCollData();

 private:
  MInt m_collisionModel = -1;

  LPT<nDim>* m_lpt = nullptr;

  MInt m_solverId = -1;
  MInt m_domainId = -1;
  MInt m_noDomains = -1;
  MPI_Comm m_mpiComm;

  MInt domainId() const { return m_domainId; }
  MInt solverId() const { return m_solverId; }
  MInt noDomains() const { return m_noDomains; }
  MPI_Comm mpiComm() const { return m_mpiComm; }

  std::list<collStruct<nDim>> collList;

  MInt m_noOfSubDomains[3]{};
  MFloat m_subDomainCoordOffset[3]{};
  MFloat m_subDomainSize[3]{};
  MInt m_totalSubDomains = 1;

  maia::lpt::subDomainCollector<nDim>* subDomainContent;
  maia::lpt::subDomainCollectorEllipsoid<nDim>* subDomainContentEllipsoid;

  MInt m_outputStep = 50;
  MFloat m_offset = -1000.0;
  MInt m_ellipsoidCCD = 1;
  MBool m_includeEllipsoids = false;
  MFloat m_timeStep = -1;

  MFloat collisionCheck(lptParticleI, lptParticleI, MFloat);

  void collisionCheckSphereEllipsoid(lptParticleI, lptEllipticleI);

  void collisionCheckEllipsoidEllipsoid(lptEllipticleI, lptEllipticleI);

  void collisionCheckSEP(lptParticleI, lptEllipticleI);

  void collisionCheckSERetroActive(lptParticleI, lptEllipticleI);

  void collisionCheckEEProActive(lptEllipticleI, lptEllipticleI);

  void collisionCheckEERetroActive(lptEllipticleI, lptEllipticleI);

  void collisionCheckEECCD(lptEllipticleI, lptEllipticleI);

  void collisionCheckSECCD(lptParticleI, lptEllipticleI);

  MFloat checkBndryCross(MInt, lptParticleI, MFloat);

  void createMatrix_CCD(const MFloat (&MA)[16], const MFloat (&u1)[3], const MFloat (&invMB)[16], const MFloat (&u2)[3],
                        const MFloat& a_2, const MFloat& c_2, MFloat (&n)[16], MFloat (&m)[5]);

  MFloat checkStaticOverlap_CCD(const MFloat (&b_Matrix)[16], const MFloat (&A)[16], const MFloat (&m)[5],
                                const MFloat& detMinusB, MFloat& t1);

  MBool checkCSI_CCD(const MFloat (&A)[16], const MFloat (&b_Matrix)[16], const MFloat (&m)[5], const MFloat& detMinusB,
                     MFloat& t1, MFloat& t2, MFloat& positive1, MFloat& positive2, MInt& iteration);

  MFloat testU_CCD(const MFloat& u, const MFloat& t4, const MFloat& t3, const MFloat& t2, const MFloat& t1,
                   const MFloat& t0);

  MFloat checkDynamicOverlap_CCD(const MFloat (&n)[16], const MFloat (&m)[5], MFloat t1, MFloat t2);

  void writeCollEvent(MFloat collisionTime, MInt particle0, MInt particle1, MFloat diam0, MFloat diam1, MFloat dens0,
                      MFloat dens1, MFloat relVel, MFloat meanVel, MFloat* x);

  // structure for collision queue element and collision queue, used to save collision datas
  // for collective writing into a Netcdf file.

  struct collQueueElem {
    MFloat collTime;
    MFloat diam0;
    MFloat diam1;
    MFloat dens0;
    MFloat dens1;
    MFloat relVel;
    MFloat meanVel;
    MFloat collPosX;
    MFloat collPosY;
    MFloat collPosZ;
    MInt part0;
    MInt part1;
  };

  std::queue<collQueueElem> collQueue;

  struct collQueueElemEllipsoid {
    MFloat collTimeStep;
    MFloat semiMinorAxis0;
    MFloat semiMinorAxis1;
    MFloat aspectRatio0;
    MFloat aspectRatio1;
    MFloat dens0;
    MFloat dens1;
    MFloat collPosX;
    MFloat collPosY;
    MFloat collPosZ;
    MInt part0;
    MInt part1;
  };

  std::queue<collQueueElemEllipsoid> collQueueEllipsoid;
};

template <MInt nDim>
class compareParticleIds {
 public:
  explicit compareParticleIds(typename std::vector<LPTSpherical<nDim>>::const_iterator aPartId) : partId(aPartId) {}

  MBool operator()(const collStruct<nDim>& obj) {
    return ((obj.returnParticle0() == partId) || (obj.returnParticle1() == partId));
  }

 private:
  typename std::vector<LPTSpherical<nDim>>::const_iterator partId;
};

#endif // MAIA_PARTICLECOLLISION_H
