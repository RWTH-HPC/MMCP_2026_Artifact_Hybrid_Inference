// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "scratch.h"

using namespace std;

size_t Scratch::m_number_of_cells;
size_t Scratch::m_number_of_elements; // = 100;

ScratchList Scratch::m_scratchSpaces;
char* Scratch::m_totalScratch; // = new MInt[Scratch::m_number_of_elements];;
size_t Scratch::m_usedmemsize;
char* Scratch::m_nextfree;
MInt Scratch::m_object_id = 0;
char* Scratch::m_maxused;
string Scratch::m_report;

#if defined(MAIA_CLANG_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
#endif

Scratch::Scratch(MFloat size, MInt Cells) {
  m_number_of_cells = Cells;
  m_number_of_elements = floor(size * sizeof(MFloat) * m_number_of_cells) + ALIGNMENT_BOUNDARY;
  m_totalScratch = new char[m_number_of_elements];
  m_nextfree = m_totalScratch;
  m_maxused = m_totalScratch;

  m_nextfree += ((ALIGNMENT_BOUNDARY - (((uintptr_t)m_nextfree) % ALIGNMENT_BOUNDARY)) % ALIGNMENT_BOUNDARY);
  ASSERT((((uintptr_t)m_nextfree) % ALIGNMENT_BOUNDARY == 0), "Scratch memory is not aligned");

  m_report = printSelfScratch();
}


/** \brief Returns a string summing up the scratch state information.
 *
 * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
 * \date 10.05.2011, 10.06.2011, 21.10.2011
 *
 * \return a string summing up the scratch state information
 **/
string Scratch::printSelfScratch() {
  stringstream begin, end, nextfree, number, total, totalmb, used, usedmb, available, availablemb, maxmem, maxmemmb;
  const MFloat bytesInAMegaByte = 1048576.0;
  begin << (MInt*)m_totalScratch;
  end << (MInt*)getEndPointer();
  nextfree << (MInt*)m_nextfree;
  number << m_number_of_elements;
  total << getTotalMemory();
  totalmb << (MFloat)getTotalMemory() / bytesInAMegaByte;
  used << m_usedmemsize;
  usedmb << (MFloat)m_usedmemsize / bytesInAMegaByte;
  available << getAvailableMemory();
  availablemb << (MFloat)getAvailableMemory() / bytesInAMegaByte;
  maxmem << (m_maxused - m_totalScratch);
  maxmemmb << (MFloat)(m_maxused - m_totalScratch) / bytesInAMegaByte;

  string message = "\n\nScratch:";
  message += "\n-----------------------------\n";
  message += "Bytes allocated:\t\t";
  message += number.str();
  message += "\nTotal memory:\t\t\t";
  message += total.str();
  message += "\t(";
  message += totalmb.str();
  message += "MB)";
  message += "\nUsed memory:\t\t\t";
  message += used.str();
  message += "\t(";
  message += usedmb.str();
  message += "MB)";
  message += "\nFree memory:\t\t\t";
  message += available.str();
  message += "\t(";
  message += availablemb.str();
  message += "MB)";
  message += "\nMax memory:\t\t\t";
  message += maxmem.str();
  message += "\t(";
  message += maxmemmb.str();
  message += "MB)";
  message += "\nScratch start pointer:\t\t";
  message += begin.str();
  message += "\nScratch end pointer:\t\t";
  message += end.str();
  message += "\nNext free pointer:\t\t";
  message += nextfree.str();
  message += "\n\n";

  return message;
}

/** \brief Returns a string summing up the scratch state information and all scratch space elements information.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * \return a string summing up the scratch state information and all scratch space elements information
 **/
string Scratch::printSelf() {
  string message = printSelfScratch();
  ScratchList::iterator iter;

  for(iter = Scratch::m_scratchSpaces.begin(); iter != Scratch::m_scratchSpaces.end(); iter++) {
    message += (*iter)->printSelf();
    message += "\n";
  }

  message += "\n\n";
  return message;
}

/** \brief Returns a shortened string summing up the scratch space state information.
 *
 * \author Andreas Lintermann
 * \date 10.05.2011
 *
 * This function is used in the report process for tracking the occasion of maximal
 * memory usage during program execution.
 *
 * \return a shortened string summing up the scratch space state information
 **/
string Scratch::printSelfReport() {
  string message = m_report;
  message += "\n\n";
  return message;
}

#if defined(MAIA_CLANG_COMPILER)
#pragma GCC diagnostic pop
#endif
