// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "alloc.h"

#include <algorithm>
#include <iostream>
#include "UTIL/debug.h"

using namespace std;

/** \brief Deallocates all memory allocated previously by mAlloc(...)
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
void mDealloc() {
  TRACE();
  using namespace maia::alloc;

  const MInt padSize = 50;
  MInt N = 20;
  std::vector<MString> topConsumers;
  std::vector<MLong> consumption(N);
  for(MInt k = 0; k < N; k++) {
    consumption[k] = 0;
    topConsumers.emplace_back("-");
  }
  MInt count = 0;
  for(auto i = g_allocatedObjects.begin(); i != g_allocatedObjects.end(); i++) {
    MLong size = (*i)->getElementSize();
    MString name = (*i)->getName();
    auto j = i + 1;
    while((j != g_allocatedObjects.end()) && ((*j)->getObjectId() == (*i)->getObjectId())) {
      size += (*j)->getElementSize();
      j++;
      i++;
    }
    count++;
    for(MInt k = 0; k < N; k++) {
      if(size > consumption[k]) {
        for(MInt l = N - 1; l > k; l--) {
          consumption[l] = consumption[l - 1];
        }
        consumption[k] = size;
        topConsumers.insert(topConsumers.begin() + k, name);
        topConsumers.pop_back();
        k = N;
      }
    }
  }

  N = mMin(N, count);

  for(MInt k = 0; k < padSize + 24; k++) {
    m_log << "_";
  }
  m_log << endl;
  m_log << "Freeing allocated memory...";

  reverse(g_allocatedObjects.begin(), g_allocatedObjects.end());

  for(auto& g_allocatedObject : g_allocatedObjects) {
#ifdef MAIA_EXTRA_DEBUG
    m_log << "Deallocating " << g_allocatedObject->getName() << " " << std::endl;
    m_log.flush();
#endif
    delete g_allocatedObject;
  }
  g_allocatedObjects.clear();

  m_log << " finished" << endl;
  m_log << "Maximum allocated memory was " << getMemorySize(maxAllocatedBytes()) << "." << endl;
  m_log << endl << "Top " << N << " (out of " << count << ") memory consumers are:" << endl;
  MChar buffer[100 + padSize];
  MLong sum = 0;
  for(MInt k = 0; k < N; k++) {
    sum += consumption[k];
    MString name(topConsumers[k], 0, mMin(std::size_t(padSize - 2), topConsumers[k].size()));
    name += ":";
    if(padSize > (signed)name.size()) {
      name.insert(name.end(), padSize - name.size(), ' ');
    }
    MString buf0 = name;
    MString buf1 = getMemorySize(consumption[k]);
    MFloat perc = 100.0 * ((static_cast<MFloat>(consumption[k])) / (static_cast<MFloat>(maxAllocatedBytes())));
    sprintf(buffer, "%3d%s %s %s %s%6.2f%s", k + 1, ".", buf0.c_str(), buf1.c_str(), "=", perc, "%");
    m_log << buffer << endl;
  }
  MString tmp;
  tmp.insert(tmp.end(), padSize + 6, ' ');
  tmp.insert(tmp.end(), 18, '-');
  m_log << tmp << endl;
  MString name;
  if(padSize > (signed)name.size()) {
    name.insert(name.end(), padSize - name.size(), ' ');
  }
  MString buf0 = name;
  MString buf1 = getMemorySize(sum);
  MFloat perc = 100.0 * ((static_cast<MFloat>(sum)) / (static_cast<MFloat>(maxAllocatedBytes())));
  sprintf(buffer, "%s%s %s %s %s%6.2f%s", "   ", " ", buf0.c_str(), buf1.c_str(), "=", perc, "%");
  m_log << buffer << endl;

  topConsumers.clear();

  if(allocatedBytes() > 0) {
    m_log << "Uncleared memory: " << getMemorySize(allocatedBytes()) << "." << endl;
  } else {
    m_log << "All memory cleared." << endl;
  }
  for(MInt k = 0; k < padSize + 24; k++) {
    m_log << "_";
  }
  m_log << endl << endl << endl;
}

/// Return the number of allocated bytes
MLong allocatedBytes() { return maia::alloc::g_allocatedBytes; }
/// Return the maximum number of allocated bytes
MLong maxAllocatedBytes() { return maia::alloc::g_maxAllocatedBytes; }
