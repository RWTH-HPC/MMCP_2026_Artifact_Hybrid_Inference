// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef VtkIo_H
#define VtkIo_H

#include "FV/fvcartesianbndrycell.h"
#include "FV/fvcartesiancellcollector.h"
#include "FV/fvcartesiansolverxd.h"
#include "FV/fvcartesiansyseqntraits.h"
#include "GRID/cartesiangridpoint.h"
#include "GRID/cartesianpointbasedcell.h"
#include "MEMORY/collector.h"
#include "MEMORY/scratch.h"
#include "contexttypes.h"

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

template <MInt nDim, class SysEqn>
class VtkIo {
 private:
  FvCartesianSolverXD<nDim, SysEqn>* m_solver = nullptr;
  Collector<FvBndryCell<nDim, SysEqn>>* m_bndryCells = nullptr;
  Collector<PointBasedCell<nDim>>* m_extractedCells = nullptr;
  Collector<CartesianGridPoint<nDim>>* m_gridPoints = nullptr;
  List<MInt>* m_bndryCellIds = nullptr;

  MInt m_noSolver;

  static constexpr MInt m_noDirs = 2 * nDim;

 public:
  using SolverCell = FvCell;

  VtkIo(FvCartesianSolverXD<nDim, SysEqn>*);
  ~VtkIo();


  SysEqn sysEqn() const { return m_solver->m_sysEqn; };
  SysEqn& sysEqn() { return m_solver->m_sysEqn; };
  MInt solverCount() const { return m_noSolver; };
  MInt noDomains() const { return m_solver->noDomains(); }
  MInt domainId() const { return m_solver->domainId(); }
  MPI_Comm mpiComm() const { return m_solver->mpiComm(); }

  MInt writeVtuArrayParallel(MPI_File&, void*, MPI_Offset, MPI_Offset, MPI_Offset);

  void insertDataHeader(char* data, uint64_t& memsize, uint64_t& memsizeGlobal, uint64_t& offset);
  void updateDataOffsets(uint64_t memsize, uint64_t& memsizeGlobal, uint64_t& offset, uint64_t& oldMemsizeGlobal);

  MInt estimateMemorySizeSolverwise(uint64_t, ScratchSpace<uint64_t>&, uint64_t);

  MBool initializeVtkXmlOutput(const MChar*, MString, MInt, MBool, MBool);
  void writeVtkXmlOutput(const MChar* fileName);


  template <typename uint_t = uint32_t>
  void writeVtuOutputParallel(
      const MString fileName, const MString geomFileName, const MInt noSolverSpecificVars = 0,
      const MFloatScratchSpace& solverSpecificVars = MFloatScratchSpace(1, "writeVtuOutputParallel",
                                                                        "defaultParameter1"),
      const MInt noUserVars = 0,
      const MFloatScratchSpace& userVars = MFloatScratchSpace(1, "writeVtuOutputParallel", "defaultParameter2"));

  template <typename T>
  uint64_t vtuAssembleFaceStream(std::vector<T>*&, std::vector<T>*&, uint64_t&, ScratchSpace<MInt>&, MInt&,
                                 const MBool);


  /** \brief Parallel single-file VTU output using MPI I/O
   *
   *  Parameters:
   *   m_solver->m_vtuWritePointData: specifies whether the flow field is written as point or cell data
   *   m_solver->m_vtuCoordinatesThreshold: optionally specifies a bounding box to which the output domain is truncated
   *   m_solver->m_vtuLevelThreshold: optionally specifies the maximum cell level to be saved
   *
   *  Opening the .vtu file requires Paraview>=4.0, resp. VTK>=6
   *  If more than 4.2 billion points shall be written the function will rerun with 64 bit integer output
   *  For data compression compile with 'WITH_ZLIB' macro
   *
   *  For file-format reference see http://www.vtk.org/Wiki/VTK_XML_Formats
   *
   *  Do not attempt to change to data types used (e.g. float or uint64_t) to XYZ data types since this will break
   * the output!
   *
   * \author Lennart Scheiders
   * \date 01/2014
   */
};

#endif
