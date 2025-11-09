// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STG_H_
#define STG_H_

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include "enums.h"

#include "fvcartesiansolverxd.h"
#include "fvstructuredsolver.h"


template <class SolverType>
class AccessorUnstructured;
template <class SolverType>
class AccessorStructured;

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

// Tag class with some type traits
template <MInt nDim, SolverType _>
struct SolverTraits {};

template <MInt nDim>
struct SolverTraits<nDim, MAIA_FINITE_VOLUME> {
  // using Accessor = AccessorUnstructured;
  using SolverType = FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>;
};

template <MInt nDim>
struct SolverTraits<nDim, MAIA_STRUCTURED> {
  // using Accessor = AccessorStructured;
  using SolverType = FvStructuredSolver<nDim>;
};

// Accessor type trait
template <MInt nDim, SolverType SolverType_>
struct AccessorTrait {};

template <MInt nDim>
struct AccessorTrait<nDim, MAIA_FINITE_VOLUME> {
  using AccessorType = AccessorUnstructured<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>>;
};

template <MInt nDim>
struct AccessorTrait<nDim, MAIA_STRUCTURED> {
  using AccessorType = AccessorStructured<FvStructuredSolver<nDim>>;
};

/// This iterator class wrapps around std::vector<MInt>::iterator
/// Use this iterator to iterate either over m_stgId or m_bcStgId, i.e., to iterate over lists containing stg ids
template <class base_iterator, class SolverType>
class nDim_iterator_t : public base_iterator {
 public:
  using value_type = typename base_iterator::value_type;

  // Constructors (default one for creation of invalid object)
  nDim_iterator_t() = default;
  nDim_iterator_t(base_iterator it, const AccessorUnstructured<SolverType>* parent) : base_iterator(it), p(parent) {}

  // public member functions
  // deprecated, since the new convention is to forward the call to the Accessor, i.e., instead of
  // 'it.getCellId()', use 'a->getCellId(it)', when 'it' is the iterator and 'a' is a pointer to the Accessor
  value_type getCellId() const;
  value_type getStgId() const { return *(*this); }
  value_type getNghbr(MInt dir) const;

 private:
  // to access member function of enclosing object
  const AccessorUnstructured<SolverType>* p{};
};


// Base class for CRTP; Defines the interface; All functions where only declarations are given, must be
// implemented by the derived class;
// The Accessor class is created for the LES solver; It manages the solver specific accesses;

template <class Derived, class SolverType_>
class Accessor {
 private:
  //   static const SolverType solvertype = Derived::solvertyp;
  //   using SolverType = typename SolverTraits<solvertype>::SolverType;
  using SolverType = SolverType_; // typename SolverTraits<3, SolverType_>::SolverType;

  // Private constructor + friend declaration prevents derived class
  // to inherit from wrong base class and from creating directly a base object
  Accessor(SolverType_* solver, const MPI_Comm& commStg) : m_solver(solver), m_commStg(commStg) {} //= default;
  friend Derived;

 protected:
  using Storage = std::vector<MInt>;
  using iterator = std::vector<MInt>::iterator;
  using const_iterator = std::vector<MInt>::const_iterator;

 public:
  // CRTP accessor
  Derived* d() { return static_cast<Derived*>(this); }

  virtual MFloat a_coordinate(MInt, MInt) const = 0;
  virtual MInt domainId(MInt, MInt) const = 0;
  virtual MFloat& a_pvariable(MInt, MInt) = 0;
  virtual MFloat& UInfinity() const = 0;
  virtual MFloat& PInfinity() const = 0;
  virtual MFloat& rhoInfinity() const = 0;

  // Member functions
  MInt sizeBC() const { return m_noBcCells; }
  MInt sizeStg() const { return m_stgSize; }


  using nDim_citerator = nDim_iterator_t<std::vector<MInt>::const_iterator, SolverType_>;

  // Note: This is not valid c++ because the iterator-types of the derived classes are different
  //      The proper way to do this is to define a common custom iterator type and derive the specific
  //      iterator type from this common type.
  //      An other way to do this is to implement a function which returns the next element to access and
  //      nullptr at the end. For example:
  //
  //      //reset to start
  //      for(accessType* example = exampleObj->begin(); example!=nullptr; example=exampleObj->next())
  //       example->bla();
  //      }

  //

  // 'iterateSlopes' -> iterate over all cells for which slopes have been computed
  nDim_citerator iterateSlopes();
  nDim_citerator iterateSlopes_nDim_citerator_end();
  // 'iterateB1' -> iterate over first layer of boundary cells
  nDim_citerator iterateB1();
  nDim_citerator iterateB1_nDim_citerator_end();
  // 'iterateAll' -> iterate over all stg cells
  nDim_citerator iterateAll();
  nDim_citerator iterateAll_nDim_citerator_end();

  MInt getCellId(const nDim_citerator& it_) const;
  MInt getStgId(const nDim_citerator& it_) const;
  MInt getNghbr(const nDim_citerator& it_, MInt dir) const;
  MInt getNghbrStg(const nDim_citerator& it_, MInt dir) const;

  // protected: //TODO labels:FV make this protected again
 public:
  SolverType* const m_solver;
  MPI_Comm m_commStg;
  static constexpr const MInt m_nDim = 3;

  // Cell ids of all stg cells in normal container
  Storage m_stgToCellId;
  // Stg ids of boundary cells
  Storage m_bcStgId;
  // Stg ids of all stg cells (just numbers in increasing order)
  Storage m_stgId;
  // # bc cells
  MInt m_noBcCells{};
  // # stg cells; m_stgSize=m_noBcCells
  MInt m_stgSize{};
};


// Accessor for unstructured finite volume solver
template <class SolverType>
class AccessorUnstructured : public Accessor<AccessorUnstructured<SolverType>, SolverType> {
 public:
  using typename Accessor<AccessorUnstructured<SolverType>, SolverType>::nDim_citerator;
  using Accessor<AccessorUnstructured<SolverType>, SolverType>::m_stgToCellId;
  using Accessor<AccessorUnstructured<SolverType>, SolverType>::m_bcStgId;
  using Accessor<AccessorUnstructured<SolverType>, SolverType>::m_stgId;
  using Accessor<AccessorUnstructured<SolverType>, SolverType>::m_noBcCells;
  using Accessor<AccessorUnstructured<SolverType>, SolverType>::m_stgSize;
  using Accessor<AccessorUnstructured<SolverType>, SolverType>::m_solver;

  // Constructor
  explicit AccessorUnstructured(MInt* sortedBndryCellIds, MInt size, SolverType* solver, const MPI_Comm& commStg,
                                MBool cutOff)
    : Accessor<AccessorUnstructured<SolverType>, SolverType>(solver, commStg) {
    MInt myrank = 0;
    MPI_Comm_rank(commStg, &myrank);

    // First put the boundary cells into m_stgToCellId and then append the ghost cells and finally
    // the reconstruction neighbors, which haven't been appended already
    m_stgToCellId.assign(sortedBndryCellIds, sortedBndryCellIds + size);
    m_noBcCells = size;
    assert(m_stgToCellId.size() == static_cast<MUint>(m_noBcCells));
    m_stgSize = m_noBcCells;

    if(!cutOff) {
      // append ghost cells and create map from stgGhost to stgBndry
      m_stgBndry2stgGhost.resize(m_noBcCells);
      m_stgGhost2stgBndry.assign(2 * m_noBcCells, -1);
      m_stgToCellId.resize(2 * m_noBcCells);
      for(MInt id = 0; id < m_noBcCells; ++id) {
        const MInt cellId = m_stgToCellId[id];
        const MInt bndryId = solver->a_bndryId(cellId);

        if(bndryId < 0) {
          TERMM(1, "errror");
        }

        const MInt ghostCellId = m_solver->m_bndryCells->a[bndryId].m_srfcVariables[0]->m_ghostCellId;
        m_stgToCellId[m_noBcCells + id] = ghostCellId;
        m_stgBndry2stgGhost[id] = m_noBcCells + id;
        m_stgGhost2stgBndry[m_noBcCells + id] = id;
      }
    }

    m_stgSize = m_stgToCellId.size();

    mAlloc(m_nghbrMapping, m_noBcCells, m_solver->m_cells.noRecNghbrs(), "m_nghbrMapping_", AT_);

    if(!cutOff) {
      // Check if cellId has been already assigned a stgId
      using MyMap = std::map<MInt, MInt>;
      MyMap cellIdToStgId_map;
      auto cellIdToStgId = [this, &cellIdToStgId_map](MInt cellId) -> MInt {
        std::pair<MyMap::const_iterator, MBool> res = cellIdToStgId_map.insert(MyMap::value_type(cellId, m_stgSize));
        if(res.second) ++m_stgSize;
        return res.first->second;
      };
      for(MInt i = 0; i < m_stgSize; ++i) {
        cellIdToStgId_map.insert(MyMap::value_type(m_stgToCellId[i], i));
      }

      for(MInt i = 0; i < m_noBcCells; ++i) {
        for(MInt nghbr = 0; nghbr < m_solver->a_noReconstructionNeighbors(m_stgToCellId[i]); ++nghbr) {
          const MInt nghbrId = m_solver->a_reconstructionNeighborId(m_stgToCellId[i], nghbr);
          const MInt stgId = cellIdToStgId(nghbrId);

          m_nghbrMapping[i][nghbr] = stgId;

          if(m_stgToCellId.size() <= static_cast<MUint>(stgId)) {
            /* if(!cutOff){ */
            if(m_solver->a_isBndryGhostCell(nghbrId)) {
              std::cout << "NOO My rank " << myrank << "; " << nghbrId << "; " << stgId << "; "
                        << solver->a_coordinate(nghbrId, 0) << "|" << solver->a_coordinate(nghbrId, 1) << "|"
                        << solver->a_coordinate(nghbrId, 2) << std::endl;
              TERMM(1, "ERROR");
            }
            /* } */
            m_stgToCellId.push_back(nghbrId);
          }
        }
      }

      if(m_stgToCellId.size() != static_cast<MUint>(m_stgSize)) {
        TERMM(1, "ERROR");
      }

      /* if(!cutOff)  */ m_stgGhost2stgBndry.resize(m_stgSize, -1);
    }

    // Create list of stg Ids of only the boundary cells, which are stored at the beginning for unstructured case
    m_bcStgId.resize(m_noBcCells);
    std::iota(std::begin(m_bcStgId), std::end(m_bcStgId), 0); // Fill with 0, 1, ..., m_noBcCells-1

    // Create list of stg Ids of all cells
    m_stgId.resize(m_stgSize);
    std::iota(std::begin(m_stgId), std::end(m_stgId), 0); // Fill with 0, 1, ..., m_stgSize-1


    // Print summary
    /* #ifndef NDEBUG */
    /*     std::cout << "--- STG INFO ---" << std::endl */
    /*               << myrank << "# bcCells/stgCells: " << m_noBcCells << " / " << m_stgSize << std::endl; */
    /* #endif */
  } // Constructor ends

  // Member functions
  MInt getNghbrMapping(MInt stgId, MInt nghbr) const { return m_nghbrMapping[stgId][nghbr]; }
  // TODO labels:FV Assumption that only one ghost cell per boundary cell!!!
  MInt getGhostIdFromStgId(MInt stgId) const {
    const MInt cellId = m_stgToCellId[stgId];
    const MInt bndryId = m_solver->a_bndryId(cellId);
    // TODO labels:FV,totest check if it is a boundary cell
    return m_solver->m_bndryCells->a[bndryId].m_srfcVariables[0]->m_ghostCellId;
  }

  /// Solver specific implementation of the interface defined in base class
  nDim_citerator iterateSlopes() { return {m_bcStgId.begin(), this}; }
  nDim_citerator iterateSlopes_nDim_citerator_end() { return {m_bcStgId.end(), this}; }

  nDim_citerator iterateB1() { return {m_bcStgId.begin(), this}; }
  nDim_citerator iterateB1_nDim_citerator_end() { return {m_bcStgId.end(), this}; }

  nDim_citerator iterateAll() { return {m_stgId.begin(), this}; }
  nDim_citerator iterateAll_nDim_citerator_end() { return {m_stgId.end(), this}; }

  /// Helper functions for the nDim_citerator as declared in base class
  MInt getCellId(const nDim_citerator& it_) const { return m_stgToCellId[*it_]; }
  MInt getStgId(const nDim_citerator& it_) const { return *it_; }
  MInt getNghbr(const nDim_citerator& it_, MInt dir) const {
    const MInt cellId = m_stgToCellId[*it_];
    if(m_solver->a_hasNeighbor(cellId, dir) > 0) {
      return m_solver->c_neighborId(cellId, dir);
    } else {
      return -1;
    }
  }
  //  MInt getNghbrStg(const nDim_citerator& it_, MInt dir) const //not implemented yet

  MFloat a_coordinate(MInt cellId, MInt dim) const override { return m_solver->a_coordinate(cellId, dim); }
  MInt domainId() const { return m_solver->domainId(); }
  MInt domainId(MInt, MInt) const override { TERMM(-1, "Invalid call!"); }
  MFloat& a_pvariable(MInt cellId, MInt varId) override { return m_solver->a_pvariable(cellId, varId); }
  MFloat& UInfinity() const override { return m_solver->m_UInfinity; }
  MFloat& PInfinity() const override { return m_solver->m_PInfinity; }
  MFloat& rhoInfinity() const override { return m_solver->m_rhoInfinity; }

 private:
  MInt** m_nghbrMapping = nullptr;
  // evtl. introduce this also in the unstructured accesor
  // maybe also introduce a iterator to iterate over all stgGhost cells
 public: // TODO labels:FV make it private later
  std::vector<MInt> m_stgBndry2stgGhost;
  std::vector<MInt> m_stgGhost2stgBndry;
};

// Helper class to access cellIds & stg cellIds from i,j,k indices
// 3D: x: m_nCells[2]; m_start[0]
//     y: m_nCells[1]; m_start[1]
//     z: m_nCells[0]; m_start[2]
// 2D: x: m_nCells[1]; m_start[0]
//     y: m_nCells[0]; m_start[1]
template <class SolverType>
class AccessorStructured : public Accessor<AccessorStructured<SolverType>, SolverType> {
 public:
  using Accessor<AccessorStructured<SolverType>, SolverType>::m_stgToCellId;
  using Accessor<AccessorStructured<SolverType>, SolverType>::m_bcStgId;
  using Accessor<AccessorStructured<SolverType>, SolverType>::m_stgId;
  using Accessor<AccessorStructured<SolverType>, SolverType>::m_noBcCells;
  using Accessor<AccessorStructured<SolverType>, SolverType>::m_stgSize;
  using Accessor<AccessorStructured<SolverType>, SolverType>::m_solver;
  using Accessor<AccessorStructured<SolverType>, SolverType>::m_nDim;

  //   static const SolverType solvertyp = MAIA_STRUCTURED;

  // Constructor
  // Check later if it is a good idea to pass a pointer to this function; what if the memory to which the pointer
  // points to, get deleted somewhere else!
  explicit AccessorStructured(const MInt* const start, const MInt* const end, const MInt* const nCells,
                              const MInt* const stgBoxSize, MInt noGhostLayers, FvStructuredSolver<3>* solver,
                              const MPI_Comm commStg)
    : Accessor<AccessorStructured<SolverType>, SolverType>(solver, commStg),
      m_start(start),
      m_end(end),
      m_nCells(nCells),
      m_stgBoxSize(stgBoxSize),
      m_noGhostLayers(noGhostLayers) {
    /*      if (m_nDim==2) {
            m_nCells[2] = m_nCells[1];
            m_nCells[1] = m_nCells[0];
            m_nCells[0] = 1;
          }*/

    m_stgSize = 1;
    for(MInt d = 0; d < m_nDim; ++d) {
      m_stgSize *= m_stgBoxSize[d];
    }

    m_noBcCells = 1;
    for(MInt d = 0; d < m_nDim; d++) {
      const MInt nDims = m_end[d] - m_start[d];
      m_noBcCells *= nDims;
    }

    m_stgToCellId.resize(m_stgSize);
    MInt cnt = 0;
    for(MInt k = 0; k < m_stgBoxSize[2]; ++k) {
      for(MInt j = 0; j < m_stgBoxSize[1]; ++j) {
        for(MInt i = 0; i < m_stgBoxSize[0]; ++i) {
          m_stgToCellId[cnt++] = cellIndex(i, j, k);
        }
      }
    }

    m_bcStgId.resize(m_noBcCells);
    cnt = 0;
    for(MInt k = m_start[2]; k < m_end[2]; ++k) {
      for(MInt j = m_start[1]; j < m_end[1]; ++j) {
        for(MInt i = m_start[0]; i < m_end[0]; ++i) {
          m_bcStgId[cnt++] = cellIndexBC(i, j, k);
        }
      }
    }

    m_stgId.resize(m_stgSize);
    std::iota(std::begin(m_stgId), std::end(m_stgId), 0); // Fill with 0, 1, ..., m_stgSize-1
  }                                                       // Constructor ends


  MInt cellIndex(const MInt i, const MInt j, const MInt k) const { return i + (j + k * m_nCells[1]) * m_nCells[2]; }
  MInt cellIndex(const MInt* const ijk_) const { return ijk_[0] + (ijk_[1] + ijk_[2] * m_nCells[1]) * m_nCells[2]; }
  MInt cellIndexBC(const MInt i, const MInt j, const MInt k) const {
    return i + (j + k * m_stgBoxSize[1]) * m_stgBoxSize[2];
  }
  MInt cellIndexBC(const MInt* const ijk_) const {
    return ijk_[0] + (ijk_[1] + ijk_[2] * m_stgBoxSize[1]) * m_stgBoxSize[2];
  }
  MInt start(MInt dim) const { return m_start[dim]; }
  MInt end(MInt dim) const { return m_end[dim]; }

  MInt domainId() const { return m_solver->domainId(); }
  MInt domainId(MInt, MInt) const override { TERMM(-1, "Invalid call!"); }

  /// Iterator, to iterate over a certain ijk-range; it internally keeps track of the current ijk
  /// values, while iterating over lists containing stg ids
  class nDim_citerator {
   public:
    using self_type = nDim_citerator;
    using iterator_category = std::random_access_iterator_tag;
    using value_type = MInt;
    using difference_type = MInt;
    using pointer = MInt*;
    using reference = MInt&;

    // Constructors (default one for creation of invalid object)
    nDim_citerator() = default;
    nDim_citerator(const AccessorStructured* parent, const /*pointer*/ MInt* start, const /*pointer*/ MInt* end)
      : p(parent), ijk_start(start), ijk_end(end) {
      //        std::copy_n(&start[0], m_nDim, &ijk_start[0]);
      //        std::copy_n(&end[0], m_nDim, &ijk_end[0]);
      std::copy_n(&start[0], m_nDim, &ijk[0]);
      stgId = p->cellIndexBC(ijk[0], ijk[1], ijk[2]);
    }

    // public member functions
    // deprecated, since the new convention is to forward the call to the Accessor, i.e., instead of
    // 'it.getCellId()', use 'a->getCellId(it)', when 'it' is the iterator and 'a' is a pointer to the Accessor
    value_type getCellId() const { return p->cellIndex(ijk); }
    value_type getStgId() const { return stgId; }
    value_type getijk(MInt dim) const { return ijk[dim]; }
    pointer getijk() const { return ijk; }
    // TODO labels:FV implement more efficient specializations of this function
    value_type getNghbr(MInt dir) const {
      return p->cellIndex(ijk[0] + map_[dir][0], ijk[1] + map_[dir][1], ijk[2] + map_[dir][2]);
    }
    value_type getNghbrStg(MInt dir) const {
      return p->cellIndexBC(ijk[0] + map_[dir][0], ijk[1] + map_[dir][1], ijk[2] + map_[dir][2]);
    }

    // pre increment operator, i.e. return current value and move to next pos; called for i++
    /*value_type*/ self_type operator++(MInt) {
      self_type temp = *this;
      if(ijk[0] < ijk_end[0]) {
        stgId++;
        ijk[0]++;
        return temp;
      }

      if(ijk[1] < ijk_end[1]) {
        ijk[0] = ijk_start[0];
        ijk[1]++;
        stgId = p->cellIndexBC(ijk[0], ijk[1], ijk[2]);
        return temp;
      }

      if(ijk[2] < ijk_end[2]) {
        ijk[0] = ijk_start[0];
        ijk[1] = ijk_start[1];
        ijk[2]++;
        stgId = p->cellIndexBC(ijk[0], ijk[1], ijk[2]);
        return temp;
      }
      if(flag_last) {
        return p->nDim_citerator_invalid;
      }
      flag_last = true;
      return temp;
    }

    // called for ++i
    /*value_type&*/ self_type& operator++() {
      if(ijk[0] < ijk_end[0]) {
        ijk[0]++;
        ++stgId;
        return *this;
      }

      if(ijk[1] < ijk_end[1]) {
        ijk[0] = ijk_start[0];
        ijk[1]++;
        stgId = p->cellIndexBC(ijk[0], ijk[1], ijk[2]);
        return *this;
      }

      if(ijk[2] < ijk_end[2]) {
        ijk[0] = ijk_start[0];
        ijk[1] = ijk_start[1];
        ijk[2]++;
        stgId = p->cellIndexBC(ijk[0], ijk[1], ijk[2]);
        return *this;
      }
      return const_cast<AccessorStructured*>(p)->nDim_citerator_invalid;
    }

    const MInt& operator*() { return stgId; }       // emulate as if it was a iterator
    const MInt& operator*() const { return stgId; } // emulate as if it was a iterator
    const self_type* operator->() { return this; }  // emulate as if it was a iterator
    MBool operator==(const self_type& rhs) const {
      return (stgId == *rhs && getijk(0) == rhs.getijk(0) && getijk(1) == rhs.getijk(1) && getijk(2) == rhs.getijk(2));
    }
    MBool operator!=(const self_type& rhs) const { return !(*this == rhs); }

   private:
    // to access member function of enclosing object
    const AccessorStructured* p{};
    // current stgId
    value_type stgId{};
    // current ijk indices
    pointer ijk{};
    // helper indices
    const MInt* const ijk_start = nullptr;
    const MInt* const ijk_end = nullptr;
    MBool flag_last = false;
  };

  // Helper functions for structured nDim_citerator
  nDim_citerator nDim_citerator_begin(const MInt* start, const MInt* end) const { return {this, start, end}; }
  nDim_citerator nDim_citerator_begin() const { return {this, m_start, m_end}; }
  nDim_citerator nDim_citerator_end() const { return nDim_citerator_invalid; }
  // Helper functions for structured nDim_citerator ends

  /// Solver specific implementations of the interface defined in base class
  nDim_citerator iterateSlopes() {
    MInt ii = 1;
    MInt start_[3] = {ii, m_start[1] + 1, m_start[2] + 1};
    MInt end_[3] = {ii + 1, m_end[1] - 1, m_end[2] - 1};
    return {this, &start_[0], &end_[0]};
  }

  nDim_citerator iterateSlopes_nDim_citerator_end() { return nDim_citerator_invalid; }

  nDim_citerator iterateB1() {
    MInt ii = m_noGhostLayers - 1;
    MInt start_[3] = {ii, m_start[1], m_start[2]};
    MInt end_[3] = {ii + 1, m_end[1], m_end[2]};
    return {this, &start_[0], &end_[0]};
  }

  nDim_citerator iterateB1_nDim_citerator_end() { return nDim_citerator_invalid; }

  nDim_citerator iterateAll() {
    MInt start_[3] = {0, 0, 0};
    return {this, &start_[0], m_stgBoxSize};
  }

  nDim_citerator iterateAll_nDim_citerator_end() { return nDim_citerator_invalid; }

  /// Helper functions for nDim_citerator as declared in base class
  MInt getCellId(const nDim_citerator& it_) const {
    assert(m_stgToCellId[*it_] == cellIndex(it_.getijk()));
    return m_stgToCellId[*it_];
  }
  MInt getStgId(const nDim_citerator& it_) const { return *it_; }
  // TODO labels:FV implement more efficient specializations of this function, by exploiting the fact, that neigbors
  //      in some directions differ from the current id by some fixed offset
  MInt getNghbr(const nDim_citerator& it_, MInt dir) const {
    return cellIndex(it_.getijk(0) + map_[dir][0], it_.getijk(1) + map_[dir][1], it_.getijk(2) + map_[dir][2]);
  }
  MInt getNghbrStg(const nDim_citerator& it_, MInt dir) const {
    return cellIndexBC(it_.getijk(0) + map_[dir][0], it_.getijk(1) + map_[dir][1], it_.getijk(2) + map_[dir][2]);
  }

  MFloat a_coordinate(MInt cellId, MInt dim) const override { return m_solver->m_cells->coordinates[dim][cellId]; }
  MFloat& a_pvariable(MInt cellId, MInt varId) override { return m_solver->m_cells->pvariables[varId][cellId]; }
  MFloat& UInfinity() const override { return m_solver->PV->UInfinity; }
  MFloat& PInfinity() const override { return m_solver->PV->PInfinity; }
  MFloat& rhoInfinity() const override { return m_solver->CV->rhoInfinity; }

 private:
  // starting indices of respective boundary
  const MInt* const m_start; //[3] = {0,0,0};
  // ending indices of respective boundary
  const MInt* const m_end; //[3] = {1,1,1};
  // # cells in each dimension
  const MInt* const m_nCells; //[3];
  // # stg cells in each dimension
  const MInt* const m_stgBoxSize; //[3];
  const MInt m_noGhostLayers;
  // Invalid iterator
  nDim_citerator nDim_citerator_invalid;
  //
  constexpr static const MInt map_[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
};
///////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME labels:FV
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
class MSTG;
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void saveStg(std::map<MInt, MSTG<nDim, SolverTypeR, SolverTypeL>*> stgBC,
             typename SolverTraits<nDim, SolverTypeL>::SolverType* solver);

// 1) In the constructor initialize things related to LES solver:
//    Set communicators;
//    Initialize iterator of LES solver;
//    Determine size of m_stgVariables (it must include space for gradient determination)
//    and allocate m_stgVariables;
// template paramters are the solver type of the RANS and LES runs
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
class MSTG {
 public:
  using SolverTypeL_ = typename SolverTraits<nDim, SolverTypeL>::SolverType;

  using Accessor = typename AccessorTrait<nDim, SolverTypeL>::AccessorType;

  using self = MSTG<nDim, SolverTypeR, SolverTypeL>;
  friend void saveStg<nDim, SolverTypeR, SolverTypeL>(std::map<MInt, self*>, SolverTypeL_*);

  friend class FvBndryCndXD<nDim, FvSysEqnNS<nDim>>;

  MSTG(FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>* solver,
       MInt bcId,
       const MPI_Comm commStg,
       MInt* sortedBndryCellIds,
       MInt noStgLCells,
       MInt stgFaceNormalDir,
       MInt stgDir,
       MInt stgWallNormalDir,
       MInt wallDir,
       MBool cutOff);

  void bc7909();
  void init(MInt commStgRoot);

 protected:
  void saveStg(); // for adaptation
  void loadStg();

 private:
  // Methods
  void setSTGProperties();
  void calcReynoldsStressConvVelLengthScale();
  void generateNewEddies(MFloat*, MFloat*);
  void advanceEddies();
  void printSTGSummary();
  void calcTotalFluctuationCholesky();
  void calcEddieCoverage();
  void determinePeriodicCells();
  void setVb(MFloat*, MFloat*);

  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*> = nullptr>
  void getBoundaryLayerThickness();
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  void getBoundaryLayerThickness();
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*> = nullptr>
  void apply();
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  void apply();
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*> = nullptr>
  void extrapolateToGX(typename Accessor::nDim_citerator);
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*> = nullptr>
  void extrapolateToBoundary(typename Accessor::nDim_citerator);
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  void extrapolateToBoundary(typename Accessor::nDim_citerator) {}
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*> = nullptr>
  void calcStrainRate();
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  void calcStrainRate();
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*> = nullptr>
  void getInflowStartEnd(MFloat*, MFloat*);
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  void getInflowStartEnd(MFloat*, MFloat*);
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*> = nullptr>
  MFloat getCellSize(typename Accessor::nDim_citerator it);
  template <class _ = void, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  MFloat getCellSize(typename Accessor::nDim_citerator it);
  template <class _ = void, std::enable_if_t<SolverTypeR == MAIA_STRUCTURED, _*> = nullptr,
            std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  void readRANSProfileStg();
  template <class _ = void, std::enable_if_t<SolverTypeR == MAIA_FINITE_VOLUME, _*> = nullptr,
            std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr>
  void readRANSProfileStg();

  void loadStg(SolverTraits<nDim, MAIA_FINITE_VOLUME>*);

  SolverTypeL_* solver() const { return m_solver; }

  MFloat generate_rand();
  MFloat generate_rand_weighted();
  MFloat generate_rand_normalDist();
  MFloat get_angle(MFloat y, MFloat z);

  /* template <class _ = void, std::enable_if_t<SolverTypeR == MAIA_STRUCTURED, _*> = nullptr, */
  /*   std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr> */
  /* void resetAfterAdaptation(){} */
  /* template <class _ = void, std::enable_if_t<SolverTypeR == MAIA_FINITE_VOLUME, _*> = nullptr, */
  /*   std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*> = nullptr> */
  /* void resetAfterAdaptation(); */

  //    Accessor*;
  static constexpr const MInt m_nDim = nDim; // 3;
  SolverTypeL_* m_solver;
  MInt m_solverId;
  const MInt m_bcId; // e.g. 7909, 7910, ...

  // ToDo: generalize this to other RANS models
  typename FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>::PrimitiveVariables* PV;

  MFloat m_sutherlandConstant;
  MFloat m_sutherlandPlusOne;

  MBool m_zonal;

  MFloat* m_inflowStart = nullptr;
  MFloat* m_inflowEnd = nullptr;

  MBool m_initialRange = false;

  MPI_Comm m_commStg;
  MInt m_commStgRoot; // Check the purpose of this; by now just set to 0
  MInt m_commStgMyRank;
  MBool m_stgRootRank;
  MFloat m_stgBLT1;
  MFloat m_stgBLT2;
  MFloat m_stgBLT3;
  MFloat m_stgDelta99Inflow;
  MBool m_stgInitialStartup;
  MInt m_stgNoEddieProperties;
  MFloat** m_stgEddies = nullptr;
  MFloat* m_stgEddyStrength = nullptr;
  MInt m_stgMaxNoEddies; // to scale length of eddies
  MFloat m_stgExple;
  MFloat m_stgEddieDistribution;
  MInt m_noStgLCells; //# LES cells
  MBool m_stgLocal = false;
  MBool m_stgCreateNewEddies;
  MInt m_stgShapeFunction;
  MBool m_stgEddieLengthScales = false;
  MBool m_stgSubSup;
  MBool m_stgSupersonic;
  MInt m_stgFaceNormalDir;
  MInt m_stgWallNormalDir = 1;
  MInt m_stgDir = 0;
  MInt m_wallDir = 1;
  MInt m_periodicDir = 2;
  MFloat* m_stgLengthFactors = nullptr;
  MFloat* m_stgRSTFactors = nullptr;
  MInt m_stgMyRank;
  MFloat* m_stgVbStart = nullptr;
  MFloat* m_stgVbEnd = nullptr;
  MFloat* m_stgVbStartFS = nullptr;
  MFloat m_stgWallEnd;
  MFloat* m_stgMaxVel = nullptr;
  MFloat** m_stgLVariables = nullptr;

  // JANNIK
  MBool m_freeStreamTurbulence = false;
  MFloat m_uuFS;
  MFloat m_vvFS;
  MFloat m_wwFS;
  MFloat m_SijSijFS;
  MFloat m_BLEddieFraction;

  MBool m_isotropicTurbulence = false;
  MBool m_preliminary = false;
  MBool m_preliminaryRans2D = false;
  MBool m_newStgMethod = false;
  MFloat m_maxStreamwiseLengthscale;
  MFloat** m_stgEddieCoverage = nullptr;
  std::vector<MInt>* m_stgPeriodicCellId = nullptr;
  std::vector<MFloat> m_stgWallNormalLocations;
  std::vector<MFloat> m_stgGlobalWallNormalLocations;
  MInt m_stgGlobalNoWallNormalLocations = F0;
  MInt* m_stgGlobalNoPeriodicLocations = nullptr;
  MInt* m_stgPeriodicIndex = nullptr;

  MBool m_cylinderTransformation = false;

  MBool m_cutOff;

  // Limiters
  const MFloat epss = 1e-34;
  const MFloat eps = 1e-16;
  const MFloat epsl = 1e-13;

  // Manual parameters
  const MFloat c_mu = 0.09;
  const MFloat a1 = 1 / sqrt(c_mu);

  const MFloat timsm = 0.3; // Smoothing factor for new eddies (in time)
  const MFloat aniso = 1.0; // Anisotropic, clustered eddies

  Accessor* a;
  Accessor* LES = a;

 public:
  struct STG;

  // random number for eddy position
  std::mt19937_64 m_PRNGEddy;
  MInt m_randomEddySeed;
  inline std::mt19937_64& randomEddyPosGenerator() { return m_PRNGEddy; }
};


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
struct MSTG<nDim, SolverTypeR, SolverTypeL>::STG {
  static constexpr const MInt AVG_U = 0;
  static constexpr const MInt AVG_V = 1;
  static constexpr const MInt AVG_W = 2;
  static constexpr const MInt AVG_UU[3] = {0, 1, 2};
  static constexpr const MInt AVG_RHO = nDim;
  static constexpr const MInt AVG_P = nDim + 1;
  static constexpr const MInt FLUC_U = nDim + 2;
  static constexpr const MInt FLUC_V = nDim + 3;
  static constexpr const MInt FLUC_W = nDim + 4;
  static constexpr const MInt NU_T = nDim + 5;
  static constexpr const MInt SIJSIJ = nDim + 6;
  static constexpr const MInt LENGTH_SCALE = nDim + 7;
  static constexpr const MInt FLUC_UU = nDim + 8;
  static constexpr const MInt FLUC_UV = nDim + 9;
  static constexpr const MInt FLUC_UW = nDim + 10;
  static constexpr const MInt FLUC_VV = nDim + 11;
  static constexpr const MInt FLUC_VW = nDim + 12;
  static constexpr const MInt FLUC_WW = nDim + 13;
  static constexpr const MInt LENGTH_X = nDim + 14;
  static constexpr const MInt LENGTH_Y = nDim + 15;
  static constexpr const MInt LENGTH_Z = nDim + 16;
  static constexpr const MInt S11 = nDim + 17;
  static constexpr const MInt S22 = nDim + 18;
  static constexpr const MInt S33 = nDim + 19;
  static constexpr const MInt S12 = nDim + 20;
  static constexpr const MInt S23 = nDim + 21;
  static constexpr const MInt S13 = nDim + 22;
  static constexpr const MInt noStgVars = nDim + 23;
};


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
constexpr const MInt MSTG<nDim, SolverTypeR, SolverTypeL>::STG::AVG_UU[];


// Standalone function to interpolate solution from structured mesh to unstructured, e.g., as an IC
// The primitive variables from the structured grid will be interpolated onto the unstructured mesh.
// Evtl. call computeConservative variables() afterwards.
// TODO: make this a member function of FvCartesianSolver
template <MInt nDim, class SysEqn>
void readRANSProfile(FvCartesianSolverXD<nDim, SysEqn>* solver) {
  TRACE();

  constexpr const MInt noVars = nDim + 2;
  MFloat zCoord = 0;
  if(nDim == 2) {
    /*! \page propertyPage1
      \section zCoordFor2DInterpolation
      <code>MInt zCoordFor2DInterpolation </code>\n
      default = <code> 0</code>\n \n
      BlaBlaBlub.\n
      Possible values are:\n
      <ul>
      <li> Float </li>
      </ul>
      Keywords: <i>INTERPOLATION, IO, FINITE_VOLUME</i>
    */
    zCoord = Context::getSolverProperty<MFloat>("zCoordFor2DInterpolation", solver->solverId(), AT_, &zCoord);
  }

  // Coordinates need to be reorderd
  //  ScratchSpace<MFloat> coords(nDim, solver->a_noCells(), AT_, "coords" );
  std::vector<std::vector<MFloat>> coords(3);
  for(MInt d = 0; d < 3; ++d) {
    coords[d].resize(solver->a_noCells());
  }

  for(MInt cellId = 0; cellId < solver->a_noCells(); ++cellId) {
    for(MInt d = 0; d < nDim; ++d) {
      coords[d][cellId] = solver->a_coordinate(cellId, d);
    }
  }
  if(nDim == 2) {
    for(MInt cellId = 0; cellId < solver->a_noCells(); ++cellId) {
      coords[2][cellId] = zCoord;
    }
  }

  // hack, since prepareInterpolation takes C-type pointer to pointer as input argument
  MInt c = 0;
  std::vector<MFloat*> coords_ptr(nDim);
  for(auto& c_ptr : coords) {
    coords_ptr[c++] = c_ptr.data();
  }

  auto* structuredInterpolation = new StructuredInterpolation<3>(solver->mpiComm());
  MInt temp[] = {solver->a_noCells(), 1, 1};
  structuredInterpolation->prepareInterpolationField(&temp[0], coords_ptr.data());

  // pvariableNames saves the conversion between the naming of the primitive variables of the FvCartesianSolver and
  // the structured solver
  std::array<MString, /*solver->PV->noVariables*/ noVars> pvariableNames;
  pvariableNames[solver->PV->U] = "u";
  pvariableNames[solver->PV->V] = "v";
  if(nDim == 3) {
    pvariableNames[solver->PV->W] = "w";
  }
  pvariableNames[solver->PV->P] = "p";
  pvariableNames[solver->PV->RHO] = "rho";
  //  ScratchSpace<MFloat> vars(nDim+2, solver->a_noCells(), AT_, "vars" );
  /*  std::vector<std::vector<MFloat>> vars(nDim+2);
    for(MInt var=0; var<solver->PV->noVariables; var++) {
      vars[var].resize(solver->a_noCells());
      structuredInterpolation->interpolateField(pvariableNames[var], &vars[var][0]);
    }

    for (MInt cellId = 0; cellId < solver->a_noCells(); ++cellId) {
      for(MInt var=0; var<solver->PV->noVariables; var++) {
        solver->a_pvariable(cellId, var) = vars[var][cellId];
      }
    }*/
  std::vector<MFloat> vars(solver->a_noCells());
  for(MInt var = 0; var < noVars /*solver->PV->noVariables*/; var++) {
    structuredInterpolation->interpolateField(pvariableNames[var], &vars[0]);
    for(MInt cellId = 0; cellId < solver->a_noCells(); ++cellId) {
      solver->a_pvariable(cellId, var) = vars[cellId];
    }
  }
  m_log << " ... FINISHED!";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// IO
/////////////////////////////////////////////////////////////////////////////////////////////////
/* Assumption in the IO is that each boundary cell has its ghost cell included as a reconstrucion
 * neighbor and thus this ghost cell appears in the stg list
 */

/*
 - The saving order of the stg cells is as follows:
   1) all non-ghost cells sorted by global id
   2) all ghost cells sorted by the global id of the respective boudnary cell
 - Save optionally the coordinates and check optionally the coordinates when loading
   the data for consistency
*/
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void saveStg(std::map<MInt, MSTG<nDim, SolverTypeR, SolverTypeL>*> stgBC,
             typename SolverTraits<nDim, SolverTypeL>::SolverType* solver) {
  /* void saveStg(std::map<MInt, MSTG<nDim, SolverTypeR, SolverTypeL>*> stgBC, SolverTypeL* solver) { */
  TRACE();

  // All ranks of the solver needs to enter this function

  //  using SolverType = typename SolverTraits<nDim, MAIA_FINITE_VOLUME>::SolverType;
  using myStg = MSTG<nDim, SolverTypeR, SolverTypeL>;
  using myMap = std::map<MInt, myStg*>;
  const MInt noVars = myStg::STG::noStgVars;

  /*! \page propertyPage1
  \section stgIOCoordinates
  <code>MInt stgIOCoordinates</code>\n
  default = <code>false</code>\n \n
  .... \n \n
0  Possible values are:
  <ul>
  <li>0,1,2</li>
  </ul>
  Keywords: <i>FINITE_VOLUME, STG</i>
  */
  MInt stgIOCoordinates = 0;
  stgIOCoordinates = Context::getSolverProperty<MInt>("stgIOCoordinates", solver->solverId(), AT_, &stgIOCoordinates);

  using namespace maia::parallel_io;
  std::stringstream filename;
  filename << solver->outputDir() << "stgRestart_" << globalTimeStep << ParallelIo::fileExt();

  m_log << "Writing restart file " << filename.str() << " ..." << std::endl;

  // The stg boundaries in stgBC are sorted by bcId, but we want to have them sorted by bndryCndId
  myMap stgBC_sorted;
  for(typename myMap::const_iterator it = stgBC.begin(); it != stgBC.end(); ++it) {
    stgBC_sorted.insert(typename myMap::value_type(it->second->m_bcId, it->second));
  }

  // Find unique list of stgBndryCnds among all ranks
  MInt noStgBCs = stgBC_sorted.size();
  std::vector<MInt> stgBCIds;
  for(typename myMap::const_iterator it = stgBC_sorted.begin(); it != stgBC_sorted.end(); ++it) {
    stgBCIds.push_back(it->first);
  }

  MPI_Allreduce(MPI_IN_PLACE, &noStgBCs, 1, MPI_INT, MPI_MAX, solver->mpiComm(), AT_, "MPI_IN_PLACE", "noStgBCs");
  std::vector<MInt> stgBCIds_global(noStgBCs * solver->noDomains());
  stgBCIds.resize(noStgBCs, -1);
  MPI_Allgather(stgBCIds.data(), noStgBCs, MPI_INT, stgBCIds_global.data(), noStgBCs, MPI_INT, solver->mpiComm(), AT_,
                "stgBCIds", "stgBCIds_global");
  std::sort(stgBCIds_global.begin(), stgBCIds_global.end());
  stgBCIds_global.erase(std::unique(stgBCIds_global.begin(), stgBCIds_global.end()), stgBCIds_global.end());
  auto itm1 = std::find(stgBCIds_global.begin(), stgBCIds_global.end(), -1);
  if(itm1 != stgBCIds_global.end()) stgBCIds_global.erase(itm1);

  // if(stgBCIds_global.size() > 0){
  // std::cout << m_solver->m_solverId << " " << solver->domainId() << " " << noStgBCs << " " << stgBCIds_global.size()
  // << std::endl;

  // MPI_Barrier(solver->mpiComm(), AT_);
  //}

  {
    ParallelIo parallelIo((filename.str()).c_str(), PIO_REPLACE, solver->mpiComm());
    parallelIo.setAttribute(static_cast<MInt>(stgBCIds_global.size()), "noStgBCs");
    parallelIo.setAttributes(&stgBCIds_global[0], "stgBCIds", stgBCIds_global.size());
  }

  // Note: It seems that it is not possible to simultanously write to the same netcdf file
  //       -> therefore ensure that all ranks involved in writing 7909 info have finished,
  //          before writing 7910 info;
  for(auto stgBCId : stgBCIds_global) {
    auto it = stgBC_sorted.find(stgBCId);
    if(it != stgBC_sorted.end()) {
      // Iterate over all stg boundaries
      //  for (typename myMap::const_iterator it = stgBC_sorted.begin(); it != stgBC_sorted.end(); ++it) {
      const MInt bcId = it->first;
      myStg* s = it->second;
      assert(bcId == s->m_bcId);
      if(bcId != s->m_bcId) // TODO labels:FV doppelt
        TERMM(1, "");

      std::stringstream stgPrefix_;
      stgPrefix_ << "stgVar" << s->m_bcId << "_";
      MString stgPrefix = stgPrefix_.str();

      if(s->m_stgLocal) {
        ParallelIo parallelIo((filename.str()).c_str(), PIO_APPEND, s->m_commStg);
        // WRITE TO FILE
        parallelIo.setAttributes(&(s->m_stgMaxNoEddies), (stgPrefix + "stgMaxNoEddies").c_str(), 1);
        parallelIo.defineArray(PIO_FLOAT, (stgPrefix + "FQeddies").c_str(),
                               s->m_stgMaxNoEddies * s->m_stgNoEddieProperties);
        parallelIo.setAttribute("FQeddies", "name", (stgPrefix + "FQeddies").c_str());
        parallelIo.setOffset(s->m_stgMaxNoEddies * s->m_stgNoEddieProperties, 0);
        parallelIo.writeArray(&(s->m_stgEddies[0][0]), (stgPrefix + "FQeddies").c_str());
      }

      /* stgIOCoordinates=0: Save internal cells and ghost cells seperately; Output is globally sorted
       *                     by globalId; globalId of ghost cell is that of its internal cell;
       *                     coordinates are not written to output
       * stgIOCoordinates=1: Same as stgIOCoordinates=0, but the coordinates of the internal and ghost
       *                     cells are also written to output, so that at reading the consistency can
       *                     be double checked;
       * stgIOCoordinates=2: Coordinates are also written to output and are used at loading to put
       *                     the data at the right place
       */

      if(stgIOCoordinates == 0 || stgIOCoordinates == 1) {
        // All ranks need to enter this if-solver

        /*
         *  if (m_stgLocal): Loop over all stg cells and save stgGlobalIds & stgGlobalIdsGhost in ordered form
         *                   Determine all window cells & halo cells
         *  Determine the position of all determined window cells in m_windowCells and send it by MPI_Alltoall
         *  Determine the position of all determined halo cells in m_haloCells and check if those halo cells are already
         * existing as window cell Send the remaining list by MPI_Alltoall; each domain now knows what data it will
         * receive in the next step Send the data to the respective domains
         */

        std::vector<std::pair<MInt /*globalId*/, MInt /*stgId*/>> stgGlobalIds;
        std::vector<std::pair<MInt /*globalId*/, MInt /*stgId*/>> stgGlobalIdsGhost;
        std::map<MInt /*cellId*/, MInt /*globalId*/> stgWindow;
        std::map<MInt /*globalId*/, MInt /*stgId*/> stgHalo;
        if(s->m_stgLocal) {
          for(typename myStg::Accessor::nDim_citerator it_ = s->a->iterateAll();
              it_ != s->a->iterateAll_nDim_citerator_end();
              ++it_) {
            const MInt cellId = s->a->getCellId(it_);
            const MInt stgId = s->a->getStgId(it_);
            if(solver->a_isHalo(cellId)) {
              const MInt globalId = solver->c_globalId(cellId);
              stgHalo.insert(std::make_pair(globalId, stgId));
            } else if(solver->a_isWindow(cellId)) {
              const MInt globalId = solver->c_globalId(cellId);
              stgWindow.insert(std::make_pair(cellId, globalId));
            }
            if(!solver->a_isBndryGhostCell(cellId)) {
              const MInt globalId = solver->c_globalId(cellId);
              stgGlobalIds.push_back(std::make_pair(globalId, stgId));
            } else /*if (solver->a_isBndryGhostCell(cellId))*/ {
              // GlobalId of ghost cell is that of its internal cell
              const MInt stgIdInternal = s->a->m_stgGhost2stgBndry[stgId];
              const MInt globalId = solver->c_globalId(s->a->m_stgToCellId[stgIdInternal]);
              stgGlobalIdsGhost.push_back(std::make_pair(globalId, stgId));
            }
          }

          // Sort the lists by the globalId
          if(std::is_sorted(stgGlobalIds.begin(), stgGlobalIds.end()))
            m_log << "Non-ghost stg cells already sorted after globalId" << std::endl;
          else {
            std::sort(stgGlobalIds.begin(), stgGlobalIds.end(),
                      [](const std::pair<MInt, MInt>& a, const std::pair<MInt, MInt>& b) { return a.first < b.first; });
          }
          if(std::is_sorted(stgGlobalIdsGhost.begin(), stgGlobalIdsGhost.end()))
            m_log << "Ghost stg cells already sorted after globalId" << std::endl;
          else {
            std::sort(stgGlobalIdsGhost.begin(), stgGlobalIdsGhost.end(),
                      [](const std::pair<MInt, MInt>& a, const std::pair<MInt, MInt>& b) { return a.first < b.first; });
          }
        } // if(s->m_stgLocal)

        // At this point all ranks should have the following data:
        //    - stgGlobalIds      --> each entry is (globalId, stgId)
        //    - stgGlobalIdsGhost --> each entry is (globalId_of_internal_cell, stgId)
        //    - stgWindow         --> each entry is map(cellId, globalId)
        //    - stgHalo           --> each entry is map(globalId, stgId)

        // Auxiliary variables needed for communication
        const MInt noNeighborDomains = solver->noNeighborDomains();
        std::vector<MInt> snghbrs(noNeighborDomains);
        for(MInt d = 0; d < noNeighborDomains; ++d) {
          snghbrs[d] = solver->neighborDomain(d);
        }

        // Determine the corresponding halo domains of the window cells
        std::vector<MInt> sendcounts(noNeighborDomains);
        std::vector<MInt> stgWindowGlobalId(stgWindow.size());
        MInt cnt = 0;
        for(MInt domain = 0; domain < noNeighborDomains; ++domain) {
          for(MInt window = 0; window < solver->noWindowCells(domain); ++window) {
            auto it_ = stgWindow.find(solver->grid().windowCell(domain, window));
            if(it_ != stgWindow.end()) {
              stgWindowGlobalId[cnt++] = it_->second;
              ++sendcounts[domain];
            }
          }
        }
        if(static_cast<MUint>(cnt) != stgWindow.size()) TERMM(-1, "Scheisse Alter!!!");

        // Communicate the window cells to the corresponding halo domains
        std::vector<MInt> globalIdsOfNghbrWindowCells =
            maia::mpi::mpiExchangePointToPoint(&stgWindowGlobalId[0], //<--Input
                                               &snghbrs[0],           //<--Input
                                               noNeighborDomains,     //<--Input
                                               &sendcounts[0],        //<--Input
                                               &snghbrs[0],           //<--Input
                                               noNeighborDomains,     //<--Input
                                               solver->mpiComm(),     //<--Input
                                               solver->domainId(),    //<--Input
                                               1);                    //<--Input

        // Erase halo cells in stgHalo map, which already exists as window cell on a different stg domain
        for(auto globalId : globalIdsOfNghbrWindowCells) {
          auto it_ = stgHalo.find(globalId);
          if(it_ != stgHalo.end()) {
            stgHalo.erase(it_);
          }
        }

        // Determine the corresponding window domains of the remaining halo cells; since the key of the stgHalo map
        // is globalId, the map should be already sorted by the globalId; store sorted globalIds in vector, such that
        // it can be sent by MPI
        if(!std::is_sorted(stgHalo.begin(), stgHalo.end(),
                           [](const std::pair<MInt, MInt>& a, const std::pair<MInt, MInt>& b) {
                             return a.first < b.first;
                           })) // comp argument to is_sorted actually not necessary
        {
          TERMM(-1, "If you read this message you are stupid!");
        }
        std::vector<MInt> haloGlobalIds(stgHalo.size());
        std::fill(sendcounts.begin(), sendcounts.end(), 0);
        MInt domain = 0;
        cnt = 0;
        for(auto& h : stgHalo) {
          haloGlobalIds[cnt++] = h.first;
          while(!(h.first < solver->domainOffset(domain + 1))) {
            ++domain;
          }
          ++sendcounts[domain]; //<-- this sendcounts is the key to the outer world
        }

        // Exchange globalIds
        std::vector<MInt> stgGlobalIdsWindow_ = maia::mpi::mpiExchangePointToPoint(&haloGlobalIds[0],
                                                                                   &snghbrs[0],
                                                                                   noNeighborDomains,
                                                                                   &sendcounts[0],
                                                                                   &snghbrs[0],
                                                                                   noNeighborDomains,
                                                                                   solver->mpiComm(),
                                                                                   solver->domainId(),
                                                                                   1);

        /////////////////////////////////////////////////////////////////////////////////////////////
        // Exchange m_stgLVariables of halo cells to their corresponding window cells
        /////////////////////////////////////////////////////////////////////////////////////////////
        // Put data to communicate (from halo cells) into send buffer tempStgVarsWS
        std::vector<MFloat> tempStgVarsWS(noVars * stgHalo.size());
        cnt = 0;
        for(auto& globalId_stgId : stgHalo) {
          const MInt stgId = globalId_stgId.second;
          for(MInt var = 0; var < myStg::STG::noStgVars; ++var) {
            tempStgVarsWS[myStg::STG::noStgVars * cnt + var] = s->m_stgLVariables[var][stgId];
          }
          ++cnt;
        }

        std::vector<MFloat> tempStgVarsWR = maia::mpi::mpiExchangePointToPoint(&tempStgVarsWS[0],
                                                                               &snghbrs[0],
                                                                               noNeighborDomains,
                                                                               &sendcounts[0],
                                                                               &snghbrs[0],
                                                                               noNeighborDomains,
                                                                               solver->mpiComm(),
                                                                               solver->domainId(),
                                                                               noVars);

        // Since stgGlobalIdsWindow has been gathered from multiple neighbor domains, they might
        // not be sorted by globalId; Using globalId as key in map will automatically sort
        std::map<MInt /*globalId*/, MInt /*stgId_inside_tempStgVarsWR*/> stgGlobalIdsWindow;
        cnt = 0;
        for(auto globalId : stgGlobalIdsWindow_) {
          stgGlobalIdsWindow.insert(std::make_pair(globalId, cnt++));
        }
        std::vector<std::pair<MInt, MInt>> stgGlobalIdsWindowV(stgGlobalIdsWindow.begin(), stgGlobalIdsWindow.end());

        /////////////////////////////////////////////////////////////////////////////////////////////
        // From here on, we will work with the following variables:
        //  tempStgVarsWR
        //  stgGlobalIdsWindowV
        //  stgGlobalIds
        //  stgGlobalIdsGhost
        /////////////////////////////////////////////////////////////////////////////////////////////

        // Save m_stgLVariables in scratch space ordered by globalId; take into account window cell variables
        // which have been sent by other domains
        MInt noStgInternalCells = stgGlobalIds.size() + stgGlobalIdsWindowV.size(); //<--
        MInt noStgGhostCells = stgGlobalIdsGhost.size();                            //<--
        std::vector<MFloat> tempStgVars(noStgInternalCells * noVars);               //<--
        std::vector<MFloat> tempStgGVars(noStgGhostCells * noVars);                 //<--
        std::vector<MInt> posInSTG(stgGlobalIds.size());
        std::vector<MInt> posInSTGWindow(stgGlobalIdsWindowV.size());
        MInt cnt1 = 0, cnt2 = 0;
        // stgGlobalIds & stgGlobalIdsWindowV should be sorted by globalId at this point
        for(MInt i = 0; i < noStgInternalCells; ++i) {
          if(stgGlobalIds[cnt1] < stgGlobalIdsWindowV[cnt2]) {
            posInSTG[cnt1++] = i;
          } else {
            posInSTGWindow[cnt2++] = i;
          }
        }
        if(static_cast<MUint>(cnt1) != stgGlobalIds.size() || static_cast<MUint>(cnt2) != stgGlobalIdsWindowV.size())
          TERMM(-1, "cnt1 != stgGlobalIds.size() or cnt2 != stgGlobalIds.size()");
        for(MInt var = 0; var < noVars; ++var) {
          for(MUint i = 0; i < posInSTG.size(); ++i) {
            const MInt stgId = stgGlobalIds[i].second;
            tempStgVars[var * noStgInternalCells + posInSTG[i]] = s->m_stgLVariables[var][stgId];
          }
          for(MUint i = 0; i < posInSTGWindow.size(); ++i) {
            const MInt stgId = stgGlobalIdsWindowV[i].second;
            tempStgVars[var * noStgInternalCells + posInSTGWindow[i]] = tempStgVarsWR[var + stgId * noVars];
          }
          cnt = 0;
          for(std::vector<std::pair<MInt, MInt>>::const_iterator it_ = stgGlobalIdsGhost.begin();
              it_ != stgGlobalIdsGhost.end();
              ++it_) {
            const MInt stgId = it_->second;
            tempStgGVars[var * noStgInternalCells + cnt] = s->m_stgLVariables[var][stgId];
            ++cnt;
          }
        }

        // Write globalIds into a linear array for IO
        std::vector<MInt> globalIds(noStgInternalCells);   //<--
        std::vector<MInt> globalIdsGhost(noStgGhostCells); //<--
        for(MUint i = 0; i < posInSTG.size(); ++i) {
          globalIds[posInSTG[i]] = stgGlobalIds[i].first;
        }
        for(MUint i = 0; i < posInSTGWindow.size(); ++i) {
          globalIds[posInSTGWindow[i]] = stgGlobalIdsWindowV[i].first;
        }
        for(MInt i = 0; i < noStgGhostCells; ++i) {
          globalIdsGhost[i] = stgGlobalIdsGhost[i].first;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // SANITY CHECK
        /////////////////////////////////////////////////////////////////////////////////////////////
        MInt noStgInternalCellsGlobal = 0;
        MInt noStgGhostCellsGlobal = 0;
        MInt* displ;
        MInt* displG;
        MInt noDomains = solver->noDomains();
        MInt* noStgInternalCellsPerRank = new MInt[noDomains];
        MInt* noStgGhostCellsPerRank = new MInt[noDomains];
        MPI_Allgather(&noStgInternalCells, 1, MPI_INT, noStgInternalCellsPerRank, 1, MPI_INT, s->m_commStg, AT_,
                      "noStgInternalCells", "noStgInternalCellsPerRank");
        MPI_Allgather(&noStgGhostCells, 1, MPI_INT, noStgGhostCellsPerRank, 1, MPI_INT, s->m_commStg, AT_,
                      "noStgGhostCells", "noStgGhostCellsPerRank");
        //      if (s->m_stgMyRank == s->m_commStgRoot) {
        displ = new MInt[noDomains];
        displG = new MInt[noDomains];
        displ[0] = 0;
        displG[0] = 0;
        for(MInt i = 1; i < noDomains; ++i) {
          displ[i] = displ[i - 1] + noStgInternalCellsPerRank[i - 1];
          displG[i] = displG[i - 1] + noStgGhostCellsPerRank[i - 1];
        }
        noStgInternalCellsGlobal = std::accumulate(noStgInternalCellsPerRank, noStgInternalCellsPerRank + noDomains, 0);
        noStgGhostCellsGlobal = std::accumulate(noStgGhostCellsPerRank, noStgGhostCellsPerRank + noDomains, 0);
        //      }

        ScratchSpace<MInt> globalIdsGlobal(noStgInternalCellsGlobal, AT_, "globalIdsGlobal");
        ScratchSpace<MInt> globalIdsGhostGlobal(noStgGhostCellsGlobal, AT_, "globalIdsGhostGlobal");
        MPI_Gatherv(&globalIds[0], noStgInternalCells, MPI_INT, &globalIdsGlobal[0], &noStgInternalCellsPerRank[0],
                    &displ[0], MPI_INT, 0, s->m_commStg, AT_, "globalIds[0]", "globalIdsGlobal[0]");
        MPI_Gatherv(&globalIdsGhost[0], noStgGhostCells, MPI_INT, &globalIdsGhostGlobal[0], &noStgGhostCellsPerRank[0],
                    &displG[0], MPI_INT, 0, s->m_commStg, AT_, "globalIdsGhost[0]", "globalIdsGhostGlobal[0]");

        // Sanity checks
        if(s->m_stgMyRank == s->m_commStgRoot) {
          /* root rank has following data:
           *    noStgInternalCellsGlobal
           *    noStgGhostCellsGlobal
           *    globalIdsGlobal
           *    globalIdsGhostGlobal
           */
          // Check if global Ids are sorted
          if(!std::is_sorted(globalIdsGlobal.begin(), globalIdsGlobal.end())) TERMM(-1, "MAIA sucks");
          if(!std::is_sorted(globalIdsGhostGlobal.begin(), globalIdsGhostGlobal.end())) TERMM(-1, "MAIA sucks");
          // Since halo cells have been skipped all global ids should be unique
          for(MInt i = 0; i < noStgInternalCellsGlobal - 1; ++i) {
            if(globalIdsGlobal[i] == globalIdsGhost[i + 1]) TERMM(-1, "MAIA really sucks");
          }
          for(MInt i = 0; i < noStgGhostCellsGlobal - 1; ++i) {
            if(globalIdsGhostGlobal[i] == globalIdsGhostGlobal[i + 1]) TERMM(-1, "MAIA really sucks");
          }
        }
        /////////////////////////////////////////////////////////////////////////////////////////////
        // WRITE TO FILE
        /////////////////////////////////////////////////////////////////////////////////////////////

        /* Following variables are written to output
         * globalIds
         * globalIdsGhost
         * tempStgVars
         * tempStgGVars
         * coordinates (optionally)
         */
        ParallelIo parallelIo((filename.str()).c_str(), PIO_APPEND, solver->mpiComm());

        // Determine data offsets for intenalStgCells and ghostStgCells
        ParallelIo::size_type stgInternalOffset, stgGhostOffset;
        ParallelIo::size_type stgGlobalNoInternalCells, stgGlobalNoGhostCells;
        parallelIo.calcOffset(noStgInternalCells, &stgInternalOffset, &stgGlobalNoInternalCells, solver->mpiComm());
        parallelIo.calcOffset(noStgGhostCells, &stgGhostOffset, &stgGlobalNoGhostCells, solver->mpiComm());

        assert(stgGlobalNoInternalCells == noStgInternalCellsGlobal);
        assert(stgGlobalNoGhostCells == noStgGhostCellsGlobal);

        // Add new variables to grid file
        parallelIo.defineArray(PIO_INT, (stgPrefix + "globalIdsInternal").c_str(), noStgInternalCells);
        parallelIo.defineArray(PIO_INT, (stgPrefix + "globalIdsGhost").c_str(), noStgGhostCells);

        // Write data to file
        parallelIo.setOffset(noStgInternalCells, stgInternalOffset);
        parallelIo.writeArray(&globalIds[0], (stgPrefix + "globalIdsInternal").c_str());
        parallelIo.setOffset(noStgGhostCells, stgGhostOffset);
        parallelIo.writeArray(&globalIdsGhost[0], (stgPrefix + "globalIdsGhost").c_str());

        // Write
        parallelIo.setOffset(noStgInternalCells, stgInternalOffset);
        for(MInt var = 0; var < myStg::STG::noStgVars; ++var) {
          std::stringstream varName;
          varName << stgPrefix + std::to_string(var);
          parallelIo.defineArray(PIO_FLOAT, (varName.str()).c_str(), noStgInternalCells);
          parallelIo.writeArray(&tempStgVars[var * noStgInternalCells], (varName.str()).c_str());
        }
        parallelIo.setOffset(noStgGhostCells, stgGhostOffset);
        for(MInt var = 0; var < myStg::STG::noStgVars; ++var) {
          std::stringstream varName;
          varName << stgPrefix + "G_" + std::to_string(var);
          parallelIo.defineArray(PIO_FLOAT, (varName.str()).c_str(), noStgGhostCells);
          parallelIo.writeArray(&tempStgGVars[var * noStgGhostCells], (varName.str()).c_str());
        }

        if(stgIOCoordinates == 1) {
          ScratchSpace<MFloat> tempStgCoord(noStgInternalCells * nDim, AT_, "tempStgCoord");
          ScratchSpace<MFloat> tempStgGCoord(noStgGhostCells * nDim, AT_, "tempStgGCoord");
          cnt = 0;
          for(auto globalId : globalIds) {
            const MInt cellId = globalId - solver->domainOffset(solver->domainId());
            for(MInt d = 0; d < nDim; ++d) {
              tempStgCoord.p[noStgInternalCells * d + cnt] = s->a->a_coordinate(cellId, d);
              ++cnt;
            }
          }
          cnt = 0;
          for(auto globalId : globalIds) {
            const MInt cellId = globalId - solver->domainOffset(solver->domainId());
            const MInt bndryId = solver->a_bndryId(cellId);
            if(bndryId < 0) TERMM(-1, "");
            const MInt ghostCellId =
                solver->a_bndryGhostCellId(bndryId,
                                           0); // solver->m_bndryCells->a[bndryId].m_srfcVariables[0]->m_ghostCellId;
            for(MInt d = 0; d < nDim; ++d) {
              tempStgGCoord.p[noStgGhostCells * d + cnt] = s->a->a_coordinate(ghostCellId, d);
              ++cnt;
            }
          }
          parallelIo.setOffset(noStgInternalCells, stgInternalOffset);
          for(MInt d = 0; d < nDim; ++d) {
            std::stringstream varName;
            varName << stgPrefix + "x" + std::to_string(d);
            parallelIo.defineArray(PIO_FLOAT, (varName.str()).c_str(), noStgInternalCells);
            parallelIo.writeArray(&tempStgCoord.p[noStgInternalCells * d], (varName.str()).c_str());
          }
          parallelIo.setOffset(noStgGhostCells, stgGhostOffset);
          for(MInt d = 0; d < nDim; ++d) {
            std::stringstream varName;
            varName << stgPrefix + "G_" + "x" + std::to_string(d);
            parallelIo.defineArray(PIO_FLOAT, (varName.str()).c_str(), noStgGhostCells);
            parallelIo.writeArray(&tempStgGCoord.p[noStgGhostCells * d], (varName.str()).c_str());
          }
        }
      } else if(stgIOCoordinates == 2) {
        // This if-solver is only executed by stg-domains

        /* Output is domainwise sorted by globalId; if all stg halo cells of current domain are
         * contained as stg window cell on some other rank, than the output should also be sorted
         * globally by the globalId; Write an attribute indicating, if output is globally sorted
         * by globalId; in case of ghostcells the globalIds are the globalIds of the respective
         * internal cells; it is also stored if cell is ghost cell; two consecutive cells in output
         * can only have the same globalId, if one of them is a ghost cell;
         */

        if(s->m_stgLocal) {
          // Create mapping commStg_rank->mpiComm_rank
          MInt comm_size;
          MPI_Comm_size(s->m_commStg, &comm_size);
          std::vector<MInt> local2GlobalRank(comm_size);
          const MInt temp = solver->domainId();
          MPI_Allgather(&temp, 1, MPI_INT, local2GlobalRank.data(), 1, MPI_INT, s->m_commStg, AT_, "solver->domainId()",
                        "local2GlobalRank");
          if(!std::is_sorted(local2GlobalRank.begin(), local2GlobalRank.end())) TERMM(1, "ERROR!");

          std::map<MInt /*cellId*/, MInt /*globalId*/> stgWindow;
          for(typename myStg::Accessor::nDim_citerator it_ = s->a->iterateAll();
              it_ != s->a->iterateAll_nDim_citerator_end();
              ++it_) {
            const MInt cellId = s->a->getCellId(it_);
            if(solver->a_isBndryGhostCell(cellId))
              continue; // TODO labels:FV this skip needs to be done evtl. also in other loops in saveSTG!!!
            if(solver->a_isWindow(cellId)) {
              const MInt globalId = solver->c_globalId(cellId);
              stgWindow.insert(std::make_pair(cellId, globalId));
            }
          }

          // Auxiliary variables needed for communication
          const MInt noNeighborDomains = solver->noNeighborDomains();
          std::vector<MInt> snghbrs; //(noNeighborDomains);
          std::vector<MInt> nghbr_commSTG2mpiComm;
          for(MInt d = 0; d < noNeighborDomains; ++d) {
            auto it_ = std::find(local2GlobalRank.begin(), local2GlobalRank.end(), solver->neighborDomain(d));
            if(it_ != local2GlobalRank.end()) {
              snghbrs.push_back(std::distance(local2GlobalRank.begin(), it_)); // solver->neighborDomain(d);
              nghbr_commSTG2mpiComm.push_back(d);
            }
          }

          // Determine the corresponding halo domainIds of the window cells
          const MInt noNeighborDomainsLocal = snghbrs.size();
          std::vector<MInt> sendcounts(noNeighborDomainsLocal);
          std::vector<MInt> stgWindowGlobalId; //(stgWindow.size());
          MInt cnt = 0;
          MInt noDoubleEntries = 0;
          for(MInt domain_ = 0; domain_ < noNeighborDomainsLocal; ++domain_) {
            const MInt domain = nghbr_commSTG2mpiComm[domain_];
            for(MInt window = 0; window < solver->noWindowCells(domain); ++window) {
              auto it_ = stgWindow.find(solver->grid().windowCell(domain, window));
              if(it_ != stgWindow.end()) {
                auto it__ = std::find(stgWindowGlobalId.begin(), stgWindowGlobalId.end(), it_->second);
                if(it__ != stgWindowGlobalId.end()) {
                  noDoubleEntries++;
                }
                //              stgWindowGlobalId[cnt++] = it_->second;
                stgWindowGlobalId.push_back(it_->second);
                ++cnt;
                ++sendcounts[domain_];
              }
            }
          }

          if(static_cast<MUint>(cnt - noDoubleEntries) != stgWindow.size()) {
            std::cout << "SCHEISSE: RANK=" << solver->domainId() << ": " << cnt << "; " << noDoubleEntries << "; "
                      << stgWindow.size() << "; " << stgWindowGlobalId.size() << std::endl;
            std::cout << std::endl;
          }
          if(static_cast<MUint>(cnt) != stgWindowGlobalId.size()) TERMM(1, "NEINNNNN");

          // Communicate the window cells to the corresponding halo domains
          MInt myrank;
          MPI_Comm_rank(s->m_commStg, &myrank);

          std::vector<MInt> globalIdsOfNghbrWindowCells =
              maia::mpi::mpiExchangePointToPoint(&stgWindowGlobalId[0],  //<--Input
                                                 &snghbrs[0],            //<--Input
                                                 noNeighborDomainsLocal, //<--Input
                                                 &sendcounts[0],         //<--Input
                                                 &snghbrs[0],            //<--Input
                                                 noNeighborDomainsLocal, //<--Input
                                                 s->m_commStg,           // solver->mpiComm(),     //<--Input
                                                 myrank,                 // solver->domainId(),    //<--Input
                                                 1);                     //<--Input

          // Skip halo cells if they are contained as window cell somewhere else
          std::vector<std::pair<MInt /*globalId*/, MInt /*stgId*/>> stgGlobalIds;
          stgGlobalIds.reserve(s->a->sizeStg());
          std::vector<std::pair<MInt /*globalId*/, MInt /*isGhost*/>> stgIsGhost;
          stgIsGhost.reserve(s->a->sizeStg());
          MInt noHaloWithNoWindow = 0;
          MInt noHaloTotal = 0;
          for(typename myStg::Accessor::nDim_citerator it_ = s->a->iterateAll();
              it_ != s->a->iterateAll_nDim_citerator_end();
              ++it_) {
            const MInt cellId = s->a->getCellId(it_);
            const MInt stgId = s->a->getStgId(it_);
            /*          if (myrank==0) {
                        std::cout << "STGID=" << stgId << " : " << s->a->m_stgSize << std::endl;
                        if (stgId>10000)
                          TERMM(1, "");
                      }*/
            MInt cellId_ = cellId;
            // If cell is ghost cell, check if its internal cell is a halo cell already contained on a different domain
            if(solver->a_isBndryGhostCell(cellId)) {
              const MInt stgIdInternal = s->a->m_stgGhost2stgBndry[stgId];
              cellId_ = s->a->m_stgToCellId[stgIdInternal];
            }
            // Skip if it is a periodic cell; periodic cells are taken into account in loadStg
            if(solver->a_isPeriodic(cellId_)) continue;
            if(solver->a_isHalo(cellId_) /*&& !solver->a_isPeriodic(cellId_)*/) {
              ++noHaloTotal;
              const MInt globalId = solver->c_globalId(cellId_);
              //            if (myrank==0 && globalId==3634505)
              //              std::cout << "HALLOO YES" << stgId << std::endl;
              // TODO labels:FV,toenhance evtl. if bndry cell is halo cell, also skip its ghost cell
              // search for this globalId in globalIdsOfNghbrWindowCells; if found continue
              if(std::find(globalIdsOfNghbrWindowCells.begin(), globalIdsOfNghbrWindowCells.end(), globalId)
                 != globalIdsOfNghbrWindowCells.end()) {
                continue;
              }
              ++noHaloWithNoWindow;
              // count no cells which are not contained as window somewhere else, in that case output should be also
              // globally sorted by globaId -> check that
            }
            if(!solver->a_isBndryGhostCell(cellId)) {
              const MInt globalId = solver->c_globalId(cellId);
              //            if (myrank==0 && globalId==3634505)
              //              std::cout << "HALLOO " << stgId << std::endl;
              stgGlobalIds.push_back(std::make_pair(globalId, stgId));
              stgIsGhost.push_back(std::make_pair(globalId, 0));
            } else /*if (solver->a_isBndryGhostCell(cellId))*/ {
              // GlobalId of ghost cells are that of their internal cell
              if(stgId > 2 * s->a->m_noBcCells) {
                std::cout << "z) On rank=" << myrank << ":EOOR " << stgId << "; " << s->a->m_noBcCells << std::endl;
                TERMM(1, "ERROR");
              }
              const MInt stgIdInternal = s->a->m_stgGhost2stgBndry[stgId];
              if(stgIdInternal > s->a->m_noBcCells) {
                std::cout << "z2) On rank=" << myrank << ":EOOR " << stgIdInternal << "; " << s->a->m_noBcCells << "; "
                          << stgId << std::endl;
                TERMM(1, "ERROR");
              }
              const MInt globalId = solver->c_globalId(s->a->m_stgToCellId[stgIdInternal]);
              //            if (myrank==0 && globalId==3634505)
              //              std::cout << "HALLOO 2 " << stgId << "; " << stgIdInternal << "; " <<
              //              s->a->m_stgToCellId[stgIdInternal] <<std::endl;

              stgGlobalIds.push_back(std::make_pair(globalId, stgId));
              // Save a flag indicating if the cell is a ghost cell
              stgIsGhost.push_back(std::make_pair(globalId, 1));
            }
          }
          /*m_log*/ std::cout << "On domain " << solver->domainId() << ": " << noHaloWithNoWindow << " of "
                              << noHaloTotal << " halo cells have no window cells somewhere else" << std::endl;

          // check the following later:  if sorted domainwise, it should be also globally sorted
          // Sort the lists by the globalId
          if(std::is_sorted(stgGlobalIds.begin(), stgGlobalIds.end()))
            m_log << "stg cells already sorted after globalId" << std::endl;
          else {
            std::sort(stgIsGhost.begin(), stgIsGhost.end(),
                      [](const std::pair<MInt, MInt>& a, const std::pair<MInt, MInt>& b) { return a.first < b.first; });
            std::sort(stgGlobalIds.begin(), stgGlobalIds.end(),
                      [](const std::pair<MInt, MInt>& a, const std::pair<MInt, MInt>& b) { return a.first < b.first; });
          }

          // Save m_stgLVariables in scratch space ordered by globalId; skip halo cell already contained as window cell
          // somewhere else
          MInt noStgInternalCells = stgGlobalIds.size();                //<--
          std::vector<MFloat> tempStgVars(noStgInternalCells * noVars); //<--
          for(MInt var = 0; var < noVars; ++var) {
            cnt = 0;
            for(auto& it_ : stgGlobalIds) {
              const MInt stgId = it_.second;
              tempStgVars[var * noStgInternalCells + cnt++] = s->m_stgLVariables[var][stgId];
            }
          }

          // Write globalIds into a linear array for IO
          std::vector<MInt> globalIds(noStgInternalCells); //<--
          std::vector<MInt> isGhost(noStgInternalCells);   //<--
          cnt = 0;
          for(auto /*&*/ itGlobalId = stgGlobalIds.begin(), itIsGhost = stgIsGhost.begin();
              itGlobalId != stgGlobalIds.end() /*|| itIsGhost!=stgIsGhost.end()*/; ++itGlobalId, ++itIsGhost) {
            globalIds[cnt] = itGlobalId->first;
            isGhost[cnt] = itIsGhost->second;
            cnt++;
          }

          ///////////////////////////////////////////////////////////////////////////////////////////
          // SANITY CHECK (actually this check can be performed by only one rank)
          ///////////////////////////////////////////////////////////////////////////////////////////
          MInt* displ;
          MInt noInvolvedRanks;
          MPI_Comm_size(s->m_commStg, &noInvolvedRanks);
          MInt* noStgInternalCellsPerRank = new MInt[noInvolvedRanks];
          MPI_Allgather(&noStgInternalCells, 1, MPI_INT, noStgInternalCellsPerRank, 1, MPI_INT, s->m_commStg, AT_,
                        "noStgInternalCells", "noStgInternalCellsPerRank");
          displ = new MInt[noInvolvedRanks];
          displ[0] = 0;
          for(MInt i = 1; i < noInvolvedRanks; ++i) {
            displ[i] = displ[i - 1] + noStgInternalCellsPerRank[i - 1];
          }
          MInt noStgInternalCellsGlobal =
              std::accumulate(noStgInternalCellsPerRank, noStgInternalCellsPerRank + noInvolvedRanks, 0);

          ScratchSpace<MInt> globalIdsGlobal(noStgInternalCellsGlobal, AT_, "globalIdsGlobal");
          MPI_Allgatherv(&globalIds[0], noStgInternalCells, MPI_INT, &globalIdsGlobal[0], &noStgInternalCellsPerRank[0],
                         &displ[0], MPI_INT, s->m_commStg, AT_, "globalIds[0]", "globalIdsGlobal[0]");
          // Gather the information about halo cells without a window cell
          MInt* noHaloWithNoWindowGlobal = new MInt[noInvolvedRanks];
          MPI_Allgather(&noHaloWithNoWindow, 1, MPI_INT, &noHaloWithNoWindowGlobal[0], 1, MPI_INT, s->m_commStg, AT_,
                        "noHaloWithNoWindow", "noHaloWithNoWindowGlobal");

          // Sanity checks
          /* root rank has following data:
           *    noStgInternalCellsGlobal
           *    globalIdsGlobal
           */
          // Check if global Ids are sorted
          MInt isGloballySorted = 1;
          if(!std::is_sorted(globalIdsGlobal.begin(), globalIdsGlobal.end())) {
            isGloballySorted = 0;
            // Check if there are halo cells for which no window cells existed->only in that case it is allowed that the
            // cells are globally not sorted by globalId
            // TODO labels:FV
            //          if (std::accumulate(noHaloWithNoWindowGlobal, noHaloWithNoWindowGlobal+noInvolvedRanks, 0)==0) {
            //            TERMM(-1, "List should be globally sorted by globalId.");
            //          }
          }
          ///////////////////////////////////////////////////////////////////////////////////////////
          // WRITE TO FILE
          ///////////////////////////////////////////////////////////////////////////////////////////

          /*
           * globalIds
           * isGhost
           * tempStgVars
           * coordinates
           */
          ParallelIo parallelIo((filename.str()).c_str(), PIO_APPEND, s->m_commStg);

          parallelIo.setAttribute(isGloballySorted, (stgPrefix + "isGloballySorted").c_str());

          // Determine data offsets for intenalStgCells
          ParallelIo::size_type stgInternalOffset;
          ParallelIo::size_type stgGlobalNoInternalCells;
          parallelIo.calcOffset(noStgInternalCells, &stgInternalOffset, &stgGlobalNoInternalCells, s->m_commStg);

          assert(stgGlobalNoInternalCells == noStgInternalCellsGlobal);
          if(stgGlobalNoInternalCells != noStgInternalCellsGlobal) {
            std::cout << "ERROR " << stgGlobalNoInternalCells << "; " << noStgInternalCellsGlobal << std::endl;
            TERMM(1, "stgGlobalNoInternalCells != noStgInternalCellsGlobal");
          }

          // Add new variables to grid file
          parallelIo.defineArray(PIO_INT, (stgPrefix + "globalIdsInternal").c_str(), stgGlobalNoInternalCells);
          parallelIo.defineArray(PIO_INT, (stgPrefix + "isGhost").c_str(), stgGlobalNoInternalCells);
          //
          for(MInt var = 0; var < myStg::STG::noStgVars; ++var) {
            std::stringstream varName;
            varName << stgPrefix + std::to_string(var);
            parallelIo.defineArray(PIO_FLOAT, (varName.str()).c_str(), stgGlobalNoInternalCells);
          }
          //
          for(MInt d = 0; d < nDim; ++d) {
            std::stringstream varName;
            varName << stgPrefix + "x" + std::to_string(d);
            parallelIo.defineArray(PIO_FLOAT, (varName.str()).c_str(), stgGlobalNoInternalCells);
          }

          // Write data to file
          parallelIo.setOffset(noStgInternalCells, stgInternalOffset);
          parallelIo.writeArray(&globalIds[0], (stgPrefix + "globalIdsInternal").c_str());
          parallelIo.writeArray(&isGhost[0], (stgPrefix + "isGhost").c_str());

          // Write
          for(MInt var = 0; var < myStg::STG::noStgVars; ++var) {
            std::stringstream varName;
            varName << stgPrefix + std::to_string(var);
            parallelIo.writeArray(&tempStgVars[var * noStgInternalCells], (varName.str()).c_str());
          }

          ScratchSpace<MFloat> tempStgCoord(noStgInternalCells * nDim, AT_, "tempStgCoord");
          cnt = 0;
          for(std::vector<std::pair<MInt, MInt>>::const_iterator it_ = stgGlobalIds.begin(); it_ != stgGlobalIds.end();
              ++it_) {
            const MInt stgId = it_->second;
            const MInt cellId = s->a->m_stgToCellId[stgId]; // s->a->getCellId(stgId); // TODO labels:FV
            for(MInt d = 0; d < nDim; ++d) {
              tempStgCoord.p[noStgInternalCells * d + cnt] = s->a->a_coordinate(cellId, d);
            }

            ++cnt;
          }
          if(cnt != noStgInternalCells) TERMM(1, "errror");
          for(MInt d = 0; d < nDim; ++d) {
            std::stringstream varName;
            varName << stgPrefix + "x" + std::to_string(d);
            parallelIo.writeArray(&(tempStgCoord.p[noStgInternalCells * d]), (varName.str()).c_str());
          }
        }
      } else
        TERMM(-1, "Unknown value for property stgIOCoordinates");
    }
    MPI_Barrier(solver->mpiComm(), AT_);
    m_log << "::saveStg: Finished writing stgRestart for " << stgBCId << std::endl;
  }
}

// Load stg function is independent of RANS solver
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::loadStg(SolverTraits<nDim, MAIA_FINITE_VOLUME>*) {
  TRACE();

  m_log << "Loading stg stuff ..." << std::endl;

  /*! \page propertyPage1
    \section stgIOCoordinates
    <code>MInt stgIOCoordinates</code>\n
    default = <code>false</code>\n \n
    .... \n \n
    Possible values are:
    <ul>
    <li>0,1,2</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, STG</i>
  */
  MInt stgIOCoordinates = 0;
  stgIOCoordinates = Context::getSolverProperty<MInt>("stgIOCoordinates", m_solver->solverId(), AT_, &stgIOCoordinates);

  std::stringstream filename;
  filename << m_solver->restartDir() << "stgRestart_" << globalTimeStep << ParallelIo::fileExt();

  using namespace maia::parallel_io;

  {
    ParallelIo parallelIo((filename.str()).c_str(), PIO_READ, m_commStg /*m_solver->mpiComm()*/);

    MInt noStgBCs;
    parallelIo.getAttribute(&noStgBCs, "noStgBCs");
    std::vector<MInt> bcIds(noStgBCs);
    parallelIo.getAttribute(&bcIds[0], "stgBCIds", noStgBCs);

    // Check if current bcId is contained in stgRestart file
    if(std::find(bcIds.begin(), bcIds.end(), m_bcId) == bcIds.end())
      TERMM(-1, "Current bcId not found in stgRestart file!");
  }

  std::stringstream stgPrefix_;
  stgPrefix_ << "stgVar" << m_bcId << "_";
  MString stgPrefix = stgPrefix_.str();

  // All stg ranks read m_stgMaxNoEddies & m_stgEddies
  {
    ParallelIo parallelIo((filename.str()).c_str(), PIO_READ, m_commStg);
    MInt stgMaxNoEddies;
    parallelIo.getAttribute(&stgMaxNoEddies, (stgPrefix + "stgMaxNoEddies").c_str());
    if(stgMaxNoEddies != m_stgMaxNoEddies && !m_stgCreateNewEddies) TERMM(-1, "");
    ParallelIo::size_type stgMaxNoEddies_ = parallelIo.getArraySize((stgPrefix + "FQeddies").c_str());
    if(stgMaxNoEddies_ != stgMaxNoEddies * m_stgNoEddieProperties) TERMM(-1, "");

    parallelIo.setOffset(stgMaxNoEddies_, 0);
    parallelIo.readArray(&(m_stgEddies[0][0]), (stgPrefix + "FQeddies").c_str()); // TODO labels:FV,IO,totest check this
  }

  // if(!m_cutOff){
  if(stgIOCoordinates == 0 || stgIOCoordinates == 1) {
#if 1 // Parallel reading, i.e. not every stg rank is reading everything

    // Assumption: ghost cells are neither window nor halo cells

    // 1) All ranks read globalIdsInternal & globalIdsGhost and check for duplicates
    ParallelIo parallelIo((filename.str()).c_str(), PIO_READ, m_solver->mpiComm());
    ParallelIo::size_type noStgInternalCellsGlobal = parallelIo.getArraySize((stgPrefix + "globalIdsInternal").c_str());
    ParallelIo::size_type noStgGhostCellsGlobal = parallelIo.getArraySize((stgPrefix + "globalIdsGhost").c_str());
    std::vector<MInt> globalIdsInternalGlobal(noStgInternalCellsGlobal);
    std::vector<MInt> globalIdsGhostGlobal(noStgGhostCellsGlobal);
    parallelIo.setOffset(noStgInternalCellsGlobal, 0);
    parallelIo.readArray(globalIdsInternalGlobal.data(), (stgPrefix + "globalIdsInternal").c_str());
    parallelIo.setOffset(noStgGhostCellsGlobal, 0);
    parallelIo.readArray(globalIdsGhostGlobal.data(), (stgPrefix + "globalIdsGhost").c_str());

    if(!std::is_sorted(globalIdsInternalGlobal.begin(), globalIdsInternalGlobal.end()))
      TERMM(-1, "globalIdsInternal is not sorted in stgRestart file");
    if(!std::is_sorted(globalIdsGhostGlobal.begin(), globalIdsGhostGlobal.end()))
      TERMM(-1, "globalIdsGhost is not sorted in stgRestart file");

    // 2)
    auto lower = std::lower_bound(globalIdsInternalGlobal.begin(), globalIdsInternalGlobal.end(),
                                  m_solver->domainOffset(m_solver->domainId()));
    auto upper = std::upper_bound(globalIdsInternalGlobal.begin(), globalIdsInternalGlobal.end(),
                                  m_solver->domainOffset(m_solver->domainId() + 1));
    std::vector<MInt> offsets(2 * m_solver->noDomains());
    const MInt domOff[] = {static_cast<MInt>(std::distance(globalIdsInternalGlobal.begin(), lower)),
                           static_cast<MInt>(std::distance(globalIdsInternalGlobal.begin(), upper))};
    const MInt noStgRestartCellsLocal = domOff[1] - domOff[0];
    // Do I have any data, but I am no stg domain?
    if(!m_stgLocal && noStgRestartCellsLocal > 0) {
      m_log << "Domain " << m_solver->domainId() << " has stg data, but is not stg domain" << std::endl;
      // Data must be window data
      for(MInt i = domOff[0]; i < domOff[1]; ++i) {
        const MInt cellId = globalIdsInternalGlobal[i] - m_solver->domainOffset(m_solver->domainId());
        if(!m_solver->a_isWindow(cellId)) {
          TERMM(-1, "Domain has non-window data, even though it is no stg domain!");
        }
      }
    }
    // Sanity checks
    MPI_Gather(&domOff[0], 2, MPI_INT, offsets.data(), 2, MPI_INT, 0, m_commStg, AT_, "domOff[0]", "offsets.data()");
    if(m_commStgMyRank == m_commStgRoot) {
      if(offsets[0] != 0) TERMM(-1, "offsets[0] != 0");
      if(offsets[offsets.size() - 1] != noStgInternalCellsGlobal) TERMM(-1, "");
      for(MInt rank = 0; rank < m_solver->noDomains() - 1; ++rank) {
        if(offsets[rank * 2 + 1] != offsets[(rank + 1) * 2 + 0]) TERMM(-1, "");
      }
    }
    // Store pairs of globalId and position in restart file of the cells residing on the current rank in a map;
    // Since key is globalId and map will be sorted by key, iterating trhough this map is same as iterating through
    // globalIdsInternalGlobal
    std::map<MInt /*globalId*/, MInt /*stgIdRestart*/> globalIdStgIdRestartMap;
    for(MInt i = domOff[0]; i < domOff[1]; ++i) {
      globalIdStgIdRestartMap.insert(globalIdStgIdRestartMap.end(),
                                     std::map<MInt, MInt>::value_type(globalIdsInternalGlobal[i], i - domOff[0]));
    }

    auto lowerG = std::lower_bound(globalIdsGhostGlobal.begin(), globalIdsGhostGlobal.end(),
                                   m_solver->domainOffset(m_solver->domainId()));
    auto upperG = std::upper_bound(globalIdsGhostGlobal.begin(), globalIdsGhostGlobal.end(),
                                   m_solver->domainOffset(m_solver->domainId() + 1));
    if(std::distance(lowerG, upperG) < 1) TERMM(-1, "");
    std::vector<MInt> offsetsG(2 * m_solver->noDomains());
    const MInt domOffG[] = {static_cast<MInt>(std::distance(globalIdsGhostGlobal.begin(), lowerG)),
                            static_cast<MInt>(std::distance(globalIdsGhostGlobal.begin(), upperG))};
    MInt noStgGRestartCellsLocal = domOffG[1] - domOffG[0];
    // Sanity checks
    if(!m_stgLocal && noStgGRestartCellsLocal > 0) TERMM(-1, "");
    MPI_Gather(&domOffG[0], 2, MPI_INT, offsetsG.data(), 2, MPI_INT, 0, m_commStg, AT_, "domOffG[0]",
               "offsetsG.data()");
    if(m_commStgMyRank == m_commStgRoot) {
      if(offsetsG[0] != 0) TERMM(-1, "offsetsG[0] != 0");
      if(offsetsG[offsetsG.size() - 1] != noStgGhostCellsGlobal) TERMM(-1, "");
      for(MInt rank = 0; rank < m_solver->noDomains() - 1; ++rank) {
        if(offsetsG[rank * 2 + 1] != offsetsG[(rank + 1) * 2 + 0]) TERMM(-1, "");
      }
    }
    std::vector<MInt> globalIdsGhost(lowerG, upperG);

    //    // Store pairs of globalId and position in restart file of the cells residing on the current rank in a map
    //    std::map<MInt, MInt> globalIdGStgIdRestartMap;
    //    for (MInt i = domOffG[0]; i < domOffG[1]; ++i) {
    //      globalIdGStgIdRestartMap.insert(globalIdGStgIdRestartMap.end(), std::map<MInt,
    //      MInt>::value_type(globalIdsGhostGlobal[i] ,i-domOffG[0]));
    //    }

    // Look for all window cells on own domain
    std::map<MInt /*cellId*/, MInt /*restartStgId*/> window_cellId_restartStgId;
    //    for (MInt stgRestartIdGlobal = domOff[0]; stgRestartIdGlobal < domOff[1]; ++stgRestartIdGlobal) {
    for(auto& it : globalIdStgIdRestartMap) {
      //      const MInt cellId = globalIdsInternalGlobal[stgRestartIdGlobal] -
      //      m_solver->domainOffset(m_solver->domainId());
      const MInt cellId = it.first - m_solver->domainOffset(m_solver->domainId());
      if(m_solver->a_isWindow(cellId)) {
        window_cellId_restartStgId.insert(std::make_pair(cellId, it.second));
      }
    }

    // Auxiliary variables needed for communication
    const MInt noNeighborDomains = m_solver->noNeighborDomains();
    std::vector<MInt> snghbrs(noNeighborDomains);
    for(MInt d = 0; d < noNeighborDomains; ++d) {
      snghbrs[d] = m_solver->neighborDomain(d);
    }

    // Determine the corresponding halo domainIds of the window cells
    std::vector<MInt> sendcounts(noNeighborDomains);
    std::vector<MInt> stgWindowSortedGlobalId(window_cellId_restartStgId.size());
    std::vector<MInt> stgWindowSortedRestartStgId(window_cellId_restartStgId.size());
    MInt cnt = 0;
    for(MInt domain = 0; domain < noNeighborDomains; ++domain) {
      for(MInt window = 0; window < m_solver->noWindowCells(domain); ++window) {
        auto it = window_cellId_restartStgId.find(m_solver->grid().windowCell(domain, window));
        if(it != window_cellId_restartStgId.end()) {
          stgWindowSortedGlobalId[cnt] = m_solver->c_globalId(it->first);
          stgWindowSortedRestartStgId[cnt] = it->second;
          ++sendcounts[domain];
          ++cnt;
        }
      }
    }
    if(static_cast<MUint>(cnt) != window_cellId_restartStgId.size()) TERMM(-1, "Scheisse Alter!!!");

    // read in internal and ghost data of current rank
    ScratchSpace<MFloat> tempStgVars(noStgRestartCellsLocal, AT_, "tempStgVars");
    ScratchSpace<MFloat> tempStgGVars(noStgGRestartCellsLocal, AT_, "tempStgGVars");
    {
      //      ParallelIo parallelIo(filename, PIO_READ, m_solver->mpiComm());
      for(MInt var = 0; var < STG::noStgVars; ++var) {
        parallelIo.setOffset(noStgRestartCellsLocal, domOff[0]);
        parallelIo.readArray(tempStgVars.getPointer(), (stgPrefix + std::to_string(var)).c_str());
        parallelIo.setOffset(noStgGRestartCellsLocal, domOffG[0]);
        parallelIo.readArray(tempStgGVars.getPointer(), (stgPrefix + std::to_string(var)).c_str());
      }
    }

    // Save data to send in scratch space
    const MInt noVars = STG::noStgVars;
    MInt noWindowsToSend = window_cellId_restartStgId.size();    //<--
    std::vector<MFloat> tempStgVarsWS(noWindowsToSend * noVars); //<--
    for(MInt var = 0; var < noVars; ++var) {
      cnt = 0;
      for(auto stgIdRestart : stgWindowSortedRestartStgId) {
        tempStgVarsWS[var * noWindowsToSend + cnt++] = tempStgVars[var * noStgRestartCellsLocal + stgIdRestart];
      }
    }

    // SEND: stgWindowSortedGlobalId & tempStgVarsWS
    std::vector<MInt> globalIdsOfNghbrWindowCells =
        maia::mpi::mpiExchangePointToPoint(&stgWindowSortedGlobalId[0], //<--Input
                                           &snghbrs[0],                 //<--Input
                                           noNeighborDomains,           //<--Input
                                           &sendcounts[0],              //<--Input
                                           &snghbrs[0],                 //<--Input
                                           noNeighborDomains,           //<--Input
                                           m_solver->mpiComm(),         //<--Input
                                           m_solver->domainId(),        //<--Input
                                           1);                          //<--Input

    std::vector<MFloat> tempStgVarsWR = maia::mpi::mpiExchangePointToPoint(&tempStgVarsWS[0],
                                                                           &snghbrs[0],
                                                                           noNeighborDomains,
                                                                           &sendcounts[0],
                                                                           &snghbrs[0],
                                                                           noNeighborDomains,
                                                                           m_solver->mpiComm(),
                                                                           m_solver->domainId(),
                                                                           noVars);

    // Maps for convenience
    std::map<MInt /*globalId*/, MInt /*stgIdRestart*/>
        map2; // for halos; this is the equivalent to globalIdStgIdRestartMap for halo cells
    for(MInt i = 0; i < noStgRestartCellsLocal; ++i) {
      globalIdStgIdRestartMap.insert(std::make_pair(globalIdsInternalGlobal[i + domOff[0]], i));
    }
    for(MUint i = 0; i < globalIdsOfNghbrWindowCells.size(); ++i) {
      map2.insert(std::make_pair(globalIdsOfNghbrWindowCells[i], i));
    }

    for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
      const MInt stgId = a->getStgId(it);
      const MInt cellId = a->getCellId(it);

      if(m_solver->a_isHalo(cellId)) {
        const MInt globalId = m_solver->c_globalId(cellId);
        auto it_ = map2.find(globalId);
        if(it_ == map2.end()) TERMM(-1, "");
        const MInt stgIdRestart = it_->second;
        map2.erase(it_);
        for(MInt var = 0; var < noVars; ++var) {
          m_stgLVariables[var][stgId] = tempStgVarsWR[noVars * stgIdRestart + var];
        }
      } else if(m_solver->a_isBndryGhostCell(cellId)) {
        const MInt stgBndry = a->m_stgGhost2stgBndry[stgId];
        const MInt globalId = m_solver->c_globalId(a->m_stgToCellId[stgBndry]);
        const MInt stgIdRestart =
            std::find(globalIdsGhost.begin(), globalIdsGhost.end(), globalId) - globalIdsGhost.begin();
        for(MInt var = 0; var < noVars; ++var) {
          m_stgLVariables[var][stgId] = tempStgGVars[noVars * stgIdRestart + var];
        }
      } else {
        const MInt globalId = m_solver->c_globalId(cellId);
        auto it_ = globalIdStgIdRestartMap.find(globalId);
        if(it_ == globalIdStgIdRestartMap.end()) TERMM(-1, "");
        const MInt stgIdRestart = it_->second;
        for(MInt var = 0; var < noVars; ++var) {
          m_stgLVariables[var][stgId] = tempStgVars[noVars * stgIdRestart + var];
        }
      }
    }
    if(globalIdStgIdRestartMap.size() > 0) TERMM(-1, "");
    if(map2.size() > 0) TERMM(-1, "");

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Check for coordinates
    ///////////////////////////////////////////////////////////////////////////////////////////////
    if(m_stgLocal) {
      if(stgIOCoordinates == 1) {
        if(!parallelIo.hasDataset((stgPrefix + "x" + std::to_string(0)).c_str()))
          TERMM(-1, "Coordinates check activated, but dataset does not contain coordinates!");

        // All stg ranks read all coordinates
        //        ParallelIo parallelIo(filename, PIO_READ, m_commStg); //TODO labels:FV,IO only stg ranks should
        //        read the coordinates
        ScratchSpace<MFloat> tempStgCoord(noStgInternalCellsGlobal * nDim, AT_, "tempStgCoord");
        ScratchSpace<MFloat> tempStgGCoord(noStgGhostCellsGlobal * nDim, AT_, "tempStgGCoord");
        parallelIo.setOffset(noStgInternalCellsGlobal, 0);
        for(MInt d = 0; d < nDim; ++d) {
          std::stringstream varName;
          varName << stgPrefix + "x" + std::to_string(d);
          parallelIo.readArray(&tempStgCoord[noStgInternalCellsGlobal * d], (varName.str()).c_str());
        }
        parallelIo.setOffset(noStgGhostCellsGlobal, 0);
        for(MInt d = 0; d < nDim; ++d) {
          std::stringstream varName;
          varName << stgPrefix + "x" + std::to_string(d);
          parallelIo.readArray(&tempStgGCoord[noStgGhostCellsGlobal * d], (varName.str()).c_str());
        }

        m_log << "Building up kd-tree..." << std::endl;

        // first create kd tree for internal cells
        std::vector<Point<3>> donorPoints;
        for(MInt stgRestartId = 0; stgRestartId < noStgInternalCellsGlobal; stgRestartId++) {
          Point<3> pt(tempStgCoord[noStgInternalCellsGlobal * 0 + stgRestartId],
                      tempStgCoord[noStgInternalCellsGlobal * 1 + stgRestartId],
                      tempStgCoord[noStgInternalCellsGlobal * 2 + stgRestartId], globalIdsInternalGlobal[stgRestartId]);
          donorPoints.push_back(pt);
        }

        // build up the tree and fill it
        KDtree<nDim>* donorTreeInternal = new KDtree<nDim>(donorPoints);

        // second create kd tree for ghost cells
        donorPoints.clear();
        for(MInt stgRestartId = 0; stgRestartId < noStgGhostCellsGlobal; stgRestartId++) {
          Point<3> pt(tempStgGCoord[noStgGhostCellsGlobal * 0 + stgRestartId],
                      tempStgGCoord[noStgGhostCellsGlobal * 1 + stgRestartId],
                      tempStgGCoord[noStgGhostCellsGlobal * 2 + stgRestartId], globalIdsGhostGlobal[stgRestartId]);
          donorPoints.push_back(pt);
        }

        // build up the tree and fill it
        KDtree<nDim>* donorTreeGhost = new KDtree<nDim>(donorPoints);

        m_log << "Building up kd-tree... FINISHED!" << std::endl;

        constexpr MFloat epsDistance = 1e-4; // TODO labels:FV make it dependent on current cell size
        for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
          const MInt stgId = a->getStgId(it);
          const MInt cellId = a->getCellId(it);
          MFloat distance = 0;
          MInt globalIdFromRestart = -1;
          MInt globalId;
          Point<3> pt(a->a_coordinate(cellId, 0), a->a_coordinate(cellId, 1), a->a_coordinate(cellId, 2));

          if(m_solver->a_isBndryGhostCell(cellId)) {
            const MInt stgBndry = a->m_stgGhost2stgBndry[stgId];
            globalId = m_solver->c_globalId(a->m_stgToCellId[stgBndry]);
            globalIdFromRestart = donorTreeGhost->nearest(pt, distance);
          } else {
            globalId = m_solver->c_globalId(cellId);
            globalIdFromRestart = donorTreeInternal->nearest(pt, distance);
          }
          if(globalIdFromRestart != globalId) TERMM(-1, "");
          if(distance > epsDistance) TERMM(-1, "");
        }
        m_log << "Coordinates check consistent!" << std::endl;
      }
    }

    ////////////////////////////////////////////////////

#else // Every stg rank is reading everything
    // Assumption: ghost cells are neither window nor halo cells
    if(m_stgLocal) {
      ParallelIo parallelIo(filename, PIO_READ, m_commStg);

      // 1) All ranks read globalIdsInternal & check for duplicates
      ParallelIo::size_type noStgInternalCellsGlobal =
          parallelIo.getArraySize((stgPrefix + "globalIdsInternal").c_str());
      ParallelIo::size_type noStgGhostCellsGlobal = parallelIo.getArraySize((stgPrefix + "globalIdsGhost").c_str());
      std::vector<MInt> globalIdsInternalGlobal(noStgInternalCellsGlobal);
      std::vector<MInt> globalIdsGhostGlobal(noStgGhostCellsGlobal);
      parallelIo.setOffset(noStgInternalCellsGlobal, 0);
      parallelIo.readArray(globalIdsInternalGlobal.data(), (stgPrefix + "globalIdsInternal").c_str());
      parallelIo.setOffset(noStgGhostCellsGlobal, 0);
      parallelIo.readArray(globalIdsGhostGlobal.data(), (stgPrefix + "globalIdsGhost").c_str());

      if(!std::is_sorted(globalIdsInternalGlobal.begin(), globalIdsInternalGlobal.end()))
        TERMM(-1, "globalIdsInternal is not sorted in stgRestart file");
      if(!std::is_sorted(globalIdsGhostGlobal.begin(), globalIdsGhostGlobal.end()))
        TERMM(-1, "globalIdsGhost is not sorted in stgRestart file");

      std::map<MInt /*globalId*/, MInt /*stgRestartId*/> globalId2stgRestartId;
      std::map<MInt /*globalId*/, MInt /*stgRestartId*/> globalId2stgRestartIdG;
      for(MInt i = 0; i < noStgInternalCellsGlobal; ++i) {
        globalId2stgRestartId.insert(std::make_pair(globalIdsInternalGlobal[i], i));
      }
      for(MInt i = 0; i < noStgGhostCellsGlobal; ++i) {
        globalId2stgRestartIdG.insert(std::make_pair(globalIdsGhostGlobal[i], i));
      }

      // Create mapping stgId to stgRestartId
      // Loop over all own stg cells
      std::vector<MInt> stgId2stgRestartId.reserve(a->sizeStg());
      std::vector<MInt> stgId2stgRestartIdG;
      MInt cnt = 0, cntG = 0;
      for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
        const MInt stgId = a->getStgId(it);
        const MInt cellId = a->getCellId(it);

        if(m_solver->a_isBndryGhostCell(cellId)) {
          const MInt stgBndry = a->m_stgGhost2stgBndry[stgId];
          const MInt globalId = m_solver->c_globalId(a->m_stgToCellId[stgBndry]);
          auto it = globalId2stgRestartIdG.find(globalId);
          if(it == globalId2stgRestartIdG.end()) TERMM(-1, "");
          stgId2stgRestartIdG[cntG++] = it->second;
          globalId2stgRestartIdG.erase(it); // to avoid that it is used twice
        } else {
          const MInt globalId = m_solver->c_globalId(cellId);
          auto it = globalId2stgRestartId.find(globalId);
          if(it == globalId2stgRestartId.end()) TERMM(-1, "");
          stgId2stgRestartId[cnt++] = it->second;
          globalId2stgRestartId.erase(it); // to avoid that it is used twice
        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // Check for coordinates
      /////////////////////////////////////////////////////////////////////////////////////////////
      if(stgIOCoordinates == 1) {
        if(!parallelIo.hasDataset((stgPrefix + "x" + io_string(0)).c_str()))
          TERMM(-1, "Coordinates check activated, but dataset does not contain coordinates!");

        // All stg ranks read all coordinates
        ScratchSpace<MFloat> tempStgCoord(noStgInternalCellsGlobal * nDim, AT_, "tempStgCoord");
        ScratchSpace<MFloat> tempStgGCoord(noStgGhostCellsGlobal * nDim, AT_, "tempStgGCoord");
        parallelIo.setOffset(noStgInternalCellsGlobal, 0);
        for(MInt d = 0; d < nDim; ++d) {
          std::stringstream varName;
          varName << stgPrefix + "x" + std::to_string(d);
          parallelIo.readArray(&tempStgCoord[noStgInternalCellsGlobal * d], (varName.str()).c_str());
        }
        parallelIo.setOffset(noStgGhostCellsGlobal, 0);
        for(MInt d = 0; d < nDim; ++d) {
          std::stringstream varName;
          varName << stgPrefix + "x" + std::to_string(d);
          parallelIo.readArray(&tempStgGCoord[noStgGhostCellsGlobal * d], (varName.str()).c_str());
        }

        m_log << "Building up kd-tree..." << std::endl;

        // first create kd tree for internal cells
        std::vector<Point<3>> donorPoints;
        for(MInt stgRestartId = 0; stgRestartId < noStgInternalCellsGlobal; stgRestartId++) {
          Point<3> a(tempStgCoord[noStgInternalCellsGlobal * 0 + stgRestartId],
                     tempStgCoord[noStgInternalCellsGlobal * 1 + stgRestartId],
                     tempStgCoord[noStgInternalCellsGlobal * 2 + stgRestartId], globalIdsInternalGlobal[stgRestartId]);
          donorPoints.push_back(a);
        }

        // build up the tree and fill it
        donorTreeInternal = new KDtree<3>(donorPoints);

        // second create kd tree for ghost cells
        donorPoints.clear();
        for(MInt stgRestartId = 0; stgRestartId < noStgGhostCellsGlobal; stgRestartId++) {
          Point<3> a(tempStgGCoord[noStgGhostCellsGlobal * 0 + stgRestartId],
                     tempStgGCoord[noStgGhostCellsGlobal * 1 + stgRestartId],
                     tempStgGCoord[noStgGhostCellsGlobal * 2 + stgRestartId], globalIdsGhostGlobal[stgRestartId]);
          donorPoints.push_back(a);
        }

        // build up the tree and fill it
        donorTreeGhost = new KDtree<3>(donorPoints);

        m_log << "Building up kd-tree... FINISHED!" << std::endl;

        constexpr MFloat epsDistance = 1e-4; // TODO labels:FV make it dependent on current cell size
        cnt = 0, cntG = 0;
        for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
          const MInt stgId = a->getStgId(it);
          const MInt cellId = a->getCellId(it);

          if(m_solver->a_isBndryGhostCell(cellId)) {
            const MInt stgRestartId = stgId2stgRestartIdG[cntG++];
            MFloat distance = 0;
            for(MInt d = 0; d < nDim; ++d) {
              distance += POW2(a->a_coordinate(cellId, d) - tempStgGCoord[noStgGhostCellsGlobal * d + stgRestartId]);
            }
            if(distance > epsDistance) TERMM(-1, "");
          } else {
            const MInt stgRestartId = stgId2stgRestartIdG[cntG++];
            MFloat distance = 0;
            for(MInt d = 0; d < nDim; ++d) {
              distance += POW2(a->a_coordinate(cellId, d) - tempStgGCoord[noStgGhostCellsGlobal * d + stgRestartId]);
            }
            if(distance > epsDistance) TERMM(-1, "");
          }
        }
      } // if(stgIOCoordinates==1)

      // Read variables and populate m_stgLVariables
      ScratchSpace<MFloat> tempStgVars(noStgInternalCellsGlobal, AT_, "tempStgVars");
      ScratchSpace<MFloat> tempStgGVars(noStgGhostCellsGlobal, AT_, "tempStgGVars");

      for(MInt var = 0; var < STG::noStgVars; ++var) {
        parallelIo.setOffset(noStgInternalCellsGlobal, 0);
        parallelIo.readArray(tempStgVars.getPointer(), (stgPrefix + std::to_string(var)).c_str());
        parallelIo.setOffset(noStgGhostCellsGlobal, 0);
        parallelIo.readArray(tempStgGVars.getPointer(), (stgPrefix + std::to_string(var)).c_str());
        cnt = 0, cntG = 0;
        for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
          const MInt stgId = a->getStgId(it);
          const MInt cellId = a->getCellId(it);

          if(m_solver->a_isBndryGhostCell(cellId)) {
            m_stgLVariables[var][stgId] = tempStgGVars[stgId2stgRestartIdG[stgId]];
          } else {
            m_stgLVariables[var][stgId] = tempStgVars[stgId2stgRestartId[stgId]];
          }
        }
      }
#endif

    //////////////////////////////////////////////////////

  } else if(stgIOCoordinates == 2) {
    /* All stgLocal ranks read all data; first read list of globalIds and coordinates;
     * loop through own stg cells and find match by means of a spatial search; check if
     * corresponding globalIds match; create map of own stgId to restartStgId
     */
    if(m_stgLocal) {
      ParallelIo parallelIo((filename.str()).c_str(), PIO_READ, m_commStg);

      ParallelIo::size_type noStgInternalCellsGlobal =
          parallelIo.getArraySize((stgPrefix + "globalIdsInternal").c_str());
      std::vector<MInt> globalIdsInternalGlobal(noStgInternalCellsGlobal);
      std::vector<MInt> isGhostGlobal(noStgInternalCellsGlobal);
      parallelIo.setOffset(noStgInternalCellsGlobal, 0);
      parallelIo.readArray(globalIdsInternalGlobal.data(), (stgPrefix + "globalIdsInternal").c_str());
      parallelIo.readArray(isGhostGlobal.data(), (stgPrefix + "isGhost").c_str());

      ScratchSpace<MFloat> tempStgCoord(noStgInternalCellsGlobal * nDim, AT_, "tempStgCoord");
      for(MInt d = 0; d < nDim; ++d) {
        std::stringstream varName;
        varName << stgPrefix + "x" + std::to_string(d);
        parallelIo.readArray(&tempStgCoord[noStgInternalCellsGlobal * d], (varName.str()).c_str());
      }

      MInt myrank = 0;
      MPI_Comm_rank(m_commStg, &myrank);
#ifndef NDEBUG
      for(MInt i = 0; i < noStgInternalCellsGlobal; ++i) {
        std::cout << "Start loadStg:myrank=" << myrank << ": (" << globalIdsInternalGlobal[i] << "/" << isGhostGlobal[i]
                  << ") : " << tempStgCoord[noStgInternalCellsGlobal * 0 + i] << "|"
                  << tempStgCoord[noStgInternalCellsGlobal * 1 + i] << "|"
                  << tempStgCoord[noStgInternalCellsGlobal * 2 + i] << std::endl;
      }
#endif

      m_log << "Building up kd-tree..." << std::endl;

      // first create kd tree for internal cells
      std::vector<Point<3>> donorPoints;
      for(MInt stgRestartId = 0; stgRestartId < noStgInternalCellsGlobal; stgRestartId++) {
        Point<3> pt(tempStgCoord[noStgInternalCellsGlobal * 0 + stgRestartId],
                    tempStgCoord[noStgInternalCellsGlobal * 1 + stgRestartId],
                    tempStgCoord[noStgInternalCellsGlobal * 2 + stgRestartId], stgRestartId);
        donorPoints.push_back(pt);
      }

      // build up the tree and fill it
      auto* donorTree = new KDtree<nDim>(donorPoints);

      m_log << "Building up kd-tree... FINISHED!" << std::endl;

      // Loop over all own stg cells
      constexpr MFloat epsDistance = 1e-5; // TODO labels:FV make it dependent on current cell size
      std::vector<MInt> stgId2stgRestartId(a->sizeStg());
      for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
        const MInt stgId = a->getStgId(it);
        const MInt cellId = a->getCellId(it);
        MFloat distance = 0;
        MInt globalId;
        MInt stgRestartId;
        if(m_solver->a_isPeriodic(cellId)
           || (m_solver->a_isBndryGhostCell(cellId)
               && m_solver->a_isPeriodic(a->m_stgToCellId[a->m_stgGhost2stgBndry[stgId]]))) {
          if(m_solver->grid().periodicCartesianDir(0) > 0 || m_solver->grid().periodicCartesianDir(1) > 0)
            TERMM(1, "Periodicity only allowed in z-direction!");
          Point<3> pt(a->a_coordinate(cellId, 0),
                      a->a_coordinate(cellId, 1),
                      a->a_coordinate(cellId, 2) + m_solver->grid().periodicCartesianLength(2));
          stgRestartId = donorTree->nearest(pt, distance);
          if(distance > epsDistance) {
            Point<3> pt2(a->a_coordinate(cellId, 0),
                         a->a_coordinate(cellId, 1),
                         a->a_coordinate(cellId, 2) - m_solver->grid().periodicCartesianLength(2));
            stgRestartId = donorTree->nearest(pt2, distance);
            if(distance > epsDistance) TERMM(1, "");
          }
        } else {
          Point<3> pt(a->a_coordinate(cellId, 0), a->a_coordinate(cellId, 1), a->a_coordinate(cellId, 2));

          stgRestartId = donorTree->nearest(pt, distance);
        }

        if(m_solver->a_isBndryGhostCell(cellId)) {
          if(!isGhostGlobal[stgRestartId]) {
            TERMM(-1, "In current computation this cell is a ghost cell, but in restart file it is not!");
          }
          const MInt stgBndry = a->m_stgGhost2stgBndry[stgId];
          globalId = m_solver->c_globalId(a->m_stgToCellId[stgBndry]);
        } else {
          globalId = m_solver->c_globalId(cellId);
        }
        if(globalIdsInternalGlobal[stgRestartId] != globalId) {
          std::cout << "FAILED (" << globalId << "/" << m_solver->a_isBndryGhostCell(cellId) << "); "
                    << globalIdsInternalGlobal[stgRestartId] << "; " << m_solver->a_isBndryGhostCell(cellId) << "; "
                    << m_solver->a_isPeriodic(cellId) << "; " << a->a_coordinate(cellId, 0) << "|"
                    << a->a_coordinate(cellId, 1) << "|" << a->a_coordinate(cellId, 2) << " || "
                    << tempStgCoord[noStgInternalCellsGlobal * 0 + stgRestartId] << "|"
                    << tempStgCoord[noStgInternalCellsGlobal * 1 + stgRestartId] << "|"
                    << tempStgCoord[noStgInternalCellsGlobal * 2 + stgRestartId] << std::endl;
          TERMM(-1, "");
        }
        if(distance > epsDistance) {
          std::cout << myrank << ":FAILED2 (" << globalId << "/" << m_solver->a_isBndryGhostCell(cellId) << "); "
                    << distance << "; " << a->a_coordinate(cellId, 0) << "|" << a->a_coordinate(cellId, 1) << "|"
                    << a->a_coordinate(cellId, 2) << " || " << tempStgCoord[noStgInternalCellsGlobal * 0 + stgRestartId]
                    << "|" << tempStgCoord[noStgInternalCellsGlobal * 1 + stgRestartId] << "|"
                    << tempStgCoord[noStgInternalCellsGlobal * 2 + stgRestartId] << std::endl;
          TERMM(-1, "");
        }
        stgId2stgRestartId[stgId] = stgRestartId;
      }

      // Read variables and populate m_stgLVariables
      ScratchSpace<MFloat> tempStgVars(noStgInternalCellsGlobal, AT_, "tempStgVars");
      for(MInt var = 0; var < STG::noStgVars; ++var) {
        parallelIo.readArray(tempStgVars.getPointer(), (stgPrefix + std::to_string(var)).c_str());
        for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
          const MInt stgId = a->getStgId(it);
          m_stgLVariables[var][stgId] = tempStgVars[stgId2stgRestartId[stgId]];
        }
      }
    } // if (m_stgLocal)
  }   // else if (stgIOCoordinates==2)
      //}//!m_cutOff
}

#endif // STG_H_
