// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTDATA_H
#define POSTDATA_H

#include <algorithm>
#include "FV/fvcartesiansyseqnns.h"
#include "FV/fvcartesiansyseqnrans.h"
#include "GRID/cartesiangrid.h"
#include "cartesiansolver.h"
#include "compiler_config.h"
#include "enums.h"
#include "postcartesiancellcollector.h"

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

template <MInt nDim_>
class PostData : public maia::CartesianSolver<nDim_, PostData<nDim_>> {
 public:
  static constexpr MInt nDim = nDim_;

  using CartesianSolver = typename maia::CartesianSolver<nDim, PostData>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;
  using SolverCell = PostCell;
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  /// Collector for POST data
  maia::post::collector::PostCellCollector<nDim> m_cells;


  // used CartesianSolver
  using CartesianSolver::domainId;
  using CartesianSolver::domainOffset;
  using CartesianSolver::grid;
  using CartesianSolver::haloCellId;
  using CartesianSolver::isActive;
  using CartesianSolver::m_freeIndices;
  using CartesianSolver::m_recalcIds;
  using CartesianSolver::m_restart;
  using CartesianSolver::m_restartFile;
  using CartesianSolver::m_restartInterval;
  using CartesianSolver::m_restartOffset;
  using CartesianSolver::m_restartTimeStep;
  using CartesianSolver::m_solutionInterval;
  using CartesianSolver::m_solutionOffset;
  using CartesianSolver::m_solutionOutput;
  using CartesianSolver::m_solutionTimeSteps;
  using CartesianSolver::m_solverId;
  using CartesianSolver::maxLevel;
  using CartesianSolver::maxNoGridCells;
  using CartesianSolver::maxRefinementLevel;
  using CartesianSolver::maxUniformRefinementLevel;
  using CartesianSolver::minLevel;
  using CartesianSolver::mpiComm;
  using CartesianSolver::noHaloCells;
  using CartesianSolver::noNeighborDomains;
  using CartesianSolver::noWindowCells;
  using CartesianSolver::outputDir;
  using CartesianSolver::restartDir;
  using CartesianSolver::returnIdleRecord;
  using CartesianSolver::returnLoadRecord;
  using CartesianSolver::saveGridFlowVars;
  using CartesianSolver::solverId;
  using CartesianSolver::updateDomainInfo;
  using CartesianSolver::windowCellId;

  using Geom = Geometry<nDim>;

  PostData(MInt, GridProxy& gridProxy_, Geom& geometry_, const MPI_Comm comm);

  ~PostData(){};

  Geom* m_geometry;

  /// Access the solver's geometry
  Geom& geometry() const { return *m_geometry; }

  MInt noInternalCells() const override { return grid().noInternalCells(); }
  MFloat timeStep() const { return m_timeStep; }
  MFloat time() const override { return globalTimeStep; }

  void preTimeStep() override{};
  void postTimeStep() override{};
  MBool solutionStep() override { return true; };
  void initSolver() override;
  void finalizeInitSolver() override;
  MBool prepareRestart(MBool, MBool&) override;

  void writeRestartFile(const MBool, const MBool, const MString, MInt*) override;
  void writeRestartFile(MBool) override{};
  void saveDataFile(const MBool writeBackup,
                    const MString fileName,
                    const MInt noVars,
                    std::vector<MString>& variablesName,
                    MFloat* variables);
  void saveRestartFile(const MBool writeBackup);
  void loadRestartFile();
  void loadMeanFile(const MString fileName);
  void loadGridFlowVars(const MChar* fileName, MInt noVariables, std::vector<MString> name);
  void saveSolverSolution(const MBool, const MBool);
  void reIntAfterRestart(MBool);
  void prepareAdaptation() override;
  void postAdaptation() override;
  void finalizeAdaptation() override;
  // Sensors
  void setSensors(std::vector<std::vector<MFloat>>& sensors,
                  std::vector<MFloat>& sensorWeight,
                  std::vector<std::bitset<64>>& sensorCellFlag,
                  std::vector<MInt>& sensorSolverId) override;

  void removeChilds(const MInt) override;
  void removeCell(const MInt) override;
  void refineCell(const MInt) override;
  void swapCells(const MInt, const MInt) override;
  void swapProxy(const MInt cellId0, const MInt cellId1) override;
  MInt cellOutside(const MFloat*, const MInt, const MInt) override;
  void resizeGridMap() override;

  virtual void cleanUp(){};

  void copyGridProperties();
  void setAndAllocateSolverData(const MBool);
  void setVariableNames();

  ///\brief Returns the number of grid-cells
  MInt c_noCells() const { return grid().tree().size(); }

  MInt a_noCells() const { return m_cells.size(); }

  /// \brief Returns wether cell \p cellId has a neighbor in dir \p dir
  MInt a_hasNeighbor(const MInt cellId, const MInt dir) const {
    if(cellId > (m_cells.size() - 1)) return -1;
    return grid().tree().hasNeighbor(cellId, dir);
  }

  /// \brief Returns the coordinate of the cell from the grid().tree() \p cellId for dimension \p dir
  const MFloat& c_coordinate(const MInt cellId, const MInt dim) { return grid().tree().coordinate(cellId, dim); }

  /// \brief Returns the coordinate of the cell from the grid().tree() \p cellId for dimension \p dir
  MFloat c_coordinate(const MInt cellId, const MInt dim) const { return grid().tree().coordinate(cellId, dim); }

  /// \brief Returns the grid level of the cell \p cellId
  MInt c_level(const MInt cellId) const { return grid().tree().level(cellId); }

  // This method is necessary for the postprocessing instance to compile
  MInt a_level(const MInt) const { TERMM(1, "Make this function return something meaningful in future!"); }
  /// \brief Returns the number of children of the cell \p cellId
  MInt c_noChildren(const MInt cellId) const { return grid().tree().noChildren(cellId); }

  /// \brief Returns the grid parent id of the cell \p cellId
  MLong c_parentId(const MInt cellId) const { return grid().tree().parent(cellId); }

  /// \brief Returns the global grid id of the grid cell \p cellId
  MLong c_globalId(const MInt cellId) const { return grid().tree().globalId(cellId); }

  MLong c_childId(const MInt cellId, const MInt pos) const { return grid().tree().child(cellId, pos); }

  /// Return the number of primitive variables
  MInt noVariables() const { return m_noVariables; };

  /// \brief Returns conservative variable \p v of the cell \p cellId for variables \p varId
  MFloat& a_variable(const MInt cellId, const MInt varId) { return m_cells.variable(cellId, varId); }

  /// \brief Returns conservative variable \p v of the cell \p cellId for variables \p varId
  MFloat a_variable(const MInt cellId, const MInt varId) const { return m_cells.variable(cellId, varId); }

  MFloat& a_averagedVariable(const MInt cellId, const MInt varId) { return m_averagedVars[cellId][varId]; }

  MBool m_forceWriteRestart = false;

  MBool isMeanFile() const { return m_isMeanFile; }
  MInt fileNoVars() const { return m_fileNoVars; }
  std::vector<MString> fileVarNames() const { return m_fileVarNames; }

  /// \brief Returns IsHalo of the cell \p cellId
  MBool a_isHalo(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsHalo); }

  /// \brief Returns IsHalo of the cell \p cellId
  maia::post::cell::BitsetType::reference a_isHalo(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsHalo);
  }

  /// \brief Returns IsWindow of the cell \p cellId
  MBool a_isWindow(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsWindow); }

  /// \brief Returns IsWindow of the cell \p cellId
  maia::post::cell::BitsetType::reference a_isWindow(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsWindow);
  }

  /// \brief Returns property \p p of the cell \p cellId
  void a_resetPropertiesSolver(const MInt cellId) { m_cells.resetProperties(cellId); }

  // time related
  MFloat m_timeStep = NAN;
  MFloat m_time = NAN;

  // restart/output related properties
  MString m_outputFormat;

  std::pair<MInt, MInt> getPropertyVariableOffset(MString s) {
    MInt index = getVariablePropertyIndex(s);
    return (index > -1) ? m_variableOffset[index]
                        : std::make_pair(std::numeric_limits<MInt>::max(), std::numeric_limits<MInt>::min());
  }

  void setPropertyVariableOffset(MString s, MInt start, MInt length) {
    MInt index = getVariablePropertyIndex(s);
    if(index > -1) {
      m_variableOffset[index].first = start;
      m_variableOffset[index].second = start + length;
    }
  }

  MString getVariableName(MInt offset);

  MInt m_sourceVarsIndex;

  // Partitioning
  void setCellWeights(MFloat*) override;
  MInt noLoadTypes() const override;
  void getDefaultWeights(MFloat* weights, std::vector<MString>& names) const;
  void getLoadQuantities(MInt* const loadQuantities) const override;
  MFloat getCellLoad(const MInt cellId, const MFloat* const weights) const override;
  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) override;
  void getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings, const MBool) override;

  /// Methods to inquire solver data information
  MInt noCellDataDlb() const override { return 1; }; // communicated data: variables
  MInt cellDataTypeDlb(const MInt dataId) const override {
    if(dataId == 0) {
      return MFLOAT;
    } else {
      TERMM(1, "solverCelldataType: invalid data id " + std::to_string(dataId));
    }
  };
  MInt cellDataSizeDlb(const MInt dataId, const MInt gridCellId) override;

  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override {
    getCellDataDlb_(dataId, oldNoCells, bufferIdToCellId, data);
  };

  /// Set solver data for DLB
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override { setCellDataDlb_(dataId, data); };

 public:
  std::vector<MString> m_variablesName;

  // Optional storage location for computed/loaded mean variables
  MFloat** m_averagedVars = nullptr;

 protected:
  /// Return the number of primitive variables
  MInt getVariablePropertyIndex(MString s) {
    auto entry = std::find(m_propertyName.begin(), m_propertyName.end(), s);
    if(entry != m_propertyName.end()) {
      return std::distance(m_propertyName.begin(), entry);
    } else {
      return -1;
    }
  };

  const std::vector<MString> m_propertyName = {"primitive",
                                               "square",
                                               "kurtosis",
                                               "skewness",
                                               "statisticCombustionAnalysis",
                                               "averageVorticity",
                                               "averageSpeedOfSound",
                                               "lamb0",
                                               "du",
                                               "drho",
                                               "dp",
                                               "gradu",
                                               "ugradu",
                                               "ugradrho",
                                               "gradprho",
                                               "rhodivu",
                                               "correlation"};

  static const std::vector<std::vector<MString>> m_averageVariableName;

  std::vector<std::pair<MInt, MInt>> m_variableOffset;

  // MFloat** m_variables_ = nullptr;
  MInt m_noVariables = 0;
  MInt m_adaptationLevel = -1;
  MBool m_forceAdaptation = false;
  MBool m_adaptationSinceLastRestart = false;

  void readProperties();

  MString m_currentGridFileName{};


  // balancing:
  MBool hasSplitBalancing() const override { return true; }
  MInt m_loadBalancingReinitStage = -1;
  void resetSolver() override;
  void balancePre() override;
  void balancePost() override;
  void finalizeBalance() override;


  // Indicator if a loaded file contains computed mean variables etc.
  MBool m_isMeanFile = false;
  // Number of variables in loaded file
  MInt m_fileNoVars = -1;
  // Variable names from loaded file
  std::vector<MString> m_fileVarNames{};


 private:
  /// \brief Set the solver cell data after DLB
  template <typename dataType>
  void setCellDataDlb_(const MInt dataId, const dataType* const data) {
    TRACE();

    // Nothing to do if solver is not active
    if(!isActive()) {
      return;
    }

    // Set the variables if this is the correct reinitialization stage
    if(dataId == 0) {
      if(m_loadBalancingReinitStage == 0) {
        std::copy_n(data, noInternalCells() * noVariables(), &a_variable(0, 0));
      }
    } else {
      TERMM(1, "Unknown data id.");
    }
  }


  /// \brief Store the solver data for a given data id ordered in the given buffer for DLB
  template <typename dataType>
  void getCellDataDlb_(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                       dataType* const data) {
    TRACE();

    MInt localBufferId = 0;
    for(MInt i = 0; i < oldNoCells; i++) {
      const MInt gridCellId = bufferIdToCellId[i];

      if(gridCellId < 0) continue;

      const MInt cellId = grid().tree().grid2solver(gridCellId);
      if(cellId < 0 || cellId >= noInternalCells()) {
        continue;
      }

      switch(dataId) {
        case 0: {
          const MInt dataSize = noVariables();
          std::copy_n(&a_variable(cellId, 0), dataSize, &data[localBufferId * dataSize]);
          break;
        }
        default:
          TERMM(1, "Unknown data id.");
          break;
      }
      localBufferId++;
    }
  }

  /// \brief Static indices for accessing post variables
  /// in nDim spatial dimensions
  static constexpr std::array<MInt, nDim> getArray012() {
    IF_CONSTEXPR(nDim == 2) {
      std::array<MInt, 2> a = {0, 1};
      return a;
    }
    else {
      std::array<MInt, 3> a = {0, 1, 2};
      return a;
    }
  }
};

#endif // POSTDATA_H
