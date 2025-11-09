// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SURFACECOUPLING_H_
#define SURFACECOUPLING_H_

#include <vector>
#include "GRID/cartesiangrid.h"
#include "GRID/cartesiangridproxy.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/maiamath.h"

namespace maia::coupling {

/**
 * \brief Setting the boundary force from one surface collector to the other
 *
 * Looping over all source surfaces (e.g. cells) and apply the forces to all mapped target surfaces (e.g. bodies)
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \tparam nDim Number of dimensions
 * \tparam S    Source collector type
 * \tparam T    Target collector type
 * \tparam M    Mapping type
 * \tparam C    Conversion factor type
 *
 * \param[in] src        Source collector of surfaces
 * \param[out] tgt       Target collector of surfaces
 * \param[in] mapping    Mapping between source and target surfaces
 * \param[in] conversion Coversion factors between source and target unit system
 */
template <MInt nDim, class S, class T, class M, class C>
void setBoundaryForce(S& src, T& tgt, M& mapping, C conversion) {
  tgt.resetForce();

  for(MInt srcId = 0; srcId < src.size(); srcId++) {
    std::array<MFloat, nDim> force{};
    for(MInt n = 0; n < nDim; n++) {
      force[n] = src.force(srcId, n) * conversion.force;
    }

    for(auto&& tgtId : mapping[srcId]) {
      tgt.addForce(tgtId, force);
    }
  }
}

/**
 * \brief Setting the boundary force and torque from one surface collector to the other
 *
 * Looping over all source surfaces (e.g. cells) and apply the forces to all mapped target surfaces (e.g. bodies)
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \tparam nDim Number of dimensions
 * \tparam S    Source collector type
 * \tparam T    Target collector type
 * \tparam M    Mapping type
 * \tparam C    Conversion factor type
 *
 * \param[in] src        Source collector of surfaces
 * \param[out] tgt       Target collector of surfaces
 * \param[in] mapping    Mapping between source and target surfaces
 * \param[in] conversion Coversion factors between source and target unit system
 */
template <MInt nDim, class S, class T, class M, class C>
void setBoundaryForceAndTorque(S& src, T& tgt, M& mapping, C conversion) {
  static constexpr MInt nRot = (nDim == 3) ? 3 : 1;

  tgt.resetForce();
  tgt.resetTorque();

  for(MInt srcId = 0; srcId < src.size(); srcId++) {
    for(MInt j = 0; j < src.noForces(); j++) {
      std::array<MFloat, nDim> force{};
      for(MInt n = 0; n < nDim; n++) {
        force[n] = src.force(srcId, j, n) * conversion.force;
      }

      for(auto&& tgtId : mapping[srcId]) {
        tgt.addForce(tgtId, force);

        std::array<MFloat, nDim> r{};
        for(MInt n = 0; n < nDim; n++) {
          r[n] = (src.surfaceCenter(srcId, j, n) - tgt.a_bodyCenter(tgtId, n)) * conversion.length;
        }

        std::array<MFloat, nRot> torque{};
        if constexpr(nDim == 3) {
          torque = maia::math::cross(force, r);
        } else if constexpr(nDim == 2) {
          torque[0] = force[0] * r[1] - force[1] * r[0];
        }

        tgt.addTorque(tgtId, torque);
      }
    }
  }
}

/**
 * \brief Setting the boundary velocity from one surface collector to the other
 *
 * Looping over all source surfaces (e.g. cells) and apply the forces to all mapped target surfaces (e.g. bodies)
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \tparam nDim Number of dimensions
 * \tparam S    Source collector type
 * \tparam T    Target collector type
 * \tparam M    Mapping type
 * \tparam C    Conversion factor type
 *
 * \param[in] src        Source collector of surfaces
 * \param[out] tgt       Target collector of surfaces
 * \param[in] mapping    Mapping between source and target surfaces
 * \param[in] conversion Coversion factors between source and target unit system
 */
template <MInt nDim, class S, class T, class M, class C>
void setBoundaryVelocity(S& src, T& tgt, M& mapping, C conversion) {
  static constexpr MInt nRot = (nDim == 3) ? 3 : 1;

  if(mapping.size() == 0) {
    return;
  }

  for(MInt srcId = 0; srcId < src.size(); srcId++) {
    std::array<MFloat, nDim> velocity{};
    src.getVelocity(srcId, velocity);
    for(MInt n = 0; n < nDim; n++) {
      velocity[n] *= conversion.velocity;
    }

    std::array<MFloat, nRot> angularVelocity{};
    src.getAngularVelocity(srcId, angularVelocity);
    for(MInt n = 0; n < nRot; n++) {
      angularVelocity[n] *= conversion.angularVelocity;
    }

    for(auto&& tgtId : mapping[srcId]) {
      std::array<MFloat, nDim> r{};
      for(MInt n = 0; n < nDim; n++) {
        r[n] = (tgt.cellCenter(tgtId, n) - src.a_bodyCenter(srcId, n)) * conversion.length;
      }
      std::array<MFloat, nDim> rotationalVelocity{};
      IF_CONSTEXPR(nDim == 3) { rotationalVelocity = maia::math::cross(r, angularVelocity); }
      else IF_CONSTEXPR(nDim == 2) {
        rotationalVelocity[0] = angularVelocity[0] * r[0];
        rotationalVelocity[1] = angularVelocity[0] * r[1];
      }

      tgt.setVelocity(tgtId, velocity);
      tgt.addVelocity(tgtId, rotationalVelocity);
    }
  }
}

/**
 * \brief Simple implementation of c++20 range
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <class Iter>
class range {
  Iter m_b;
  Iter m_e;

 public:
  range(Iter b, Iter e) : m_b(b), m_e(e) {}

  Iter begin() { return m_b; }
  Iter end() { return m_e; }
};

/**
 * \brief Multi-to-multi mapping class
 *
 * This mapping class allows for multi-to-multi mappings returning a range of values for a given key
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
struct Mapping {
  std::vector<MUint> m_offsets;
  std::vector<MInt> m_mapped;

  range<std::vector<MInt>::const_iterator> get(const MUint i) const {
    MUint b = 0;
    MUint e = 0;

    if(i < m_offsets.size()) {
      b = m_offsets[i];
      e = (i == m_offsets.size() - 1) ? m_mapped.size() : m_offsets[i + 1];
    }

    std::vector<MInt>::const_iterator bb = m_mapped.begin() + b;
    std::vector<MInt>::const_iterator ee = m_mapped.begin() + e;

    return {bb, ee};
  }

  void clear() {
    m_offsets.clear();
    m_mapped.clear();
  }

  // Total count of mapped values
  std::size_t size() const { return m_mapped.size(); }

  // Insert value for key i, return value reference
  MInt& insert(const MInt i) {
    if(i < MInt(m_offsets.size() - 1)) {
      // key already exists, but is not the last one

      // if value is -1, set now!
      if(m_mapped[m_offsets[i]] == -1) {
        return m_mapped[m_offsets[i]];
      }

      for(MInt j = i + 1; j < MInt(m_offsets.size()); j++) {
        m_offsets[j]++;
      }

      m_mapped.insert(m_mapped.begin() + m_offsets[i + 1] - 1, -1);

      return m_mapped[m_offsets[i + 1] - 1];

    } else if(i == MInt(m_offsets.size() - 1)) {
      // key already exists and is the last

      m_mapped.push_back(-1);

      return m_mapped[m_mapped.size() - 1];

    } else if(i >= MInt(m_offsets.size())) {
      // key does not exist yet

      // number of keys to be added
      // none is mapping onto a value but the last one
      const MUint noNewKeys = i + 1 - m_offsets.size();

      for(MUint j = 0; j < noNewKeys; j++) {
        m_offsets.push_back(m_mapped.size());
      }

      m_mapped.push_back(-1);

      return m_mapped[m_mapped.size() - 1];
    }

    std::cerr << "Map insert shouldnt reach this case" << std::endl;
    return m_mapped[m_mapped.size() - 1];
  }

  range<std::vector<MInt>::const_iterator> operator[](const MInt i) const { return get(i); }
};


} // namespace maia::coupling
#endif
