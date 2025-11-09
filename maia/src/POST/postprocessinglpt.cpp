// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessinglpt.h"
#include <fstream>
#include <iostream>
#include "postdata.h"
#include "postprocessing.cpp"

using namespace std;

template <MInt nDim>
PostProcessingLPT<nDim>::PostProcessingLPT(MInt postprocessingId_, PostData<nDim>* data, LPT<nDim>* ppSolver_)
  : PostProcessingInterface(postprocessingId_), PostProcessing<nDim, PostProcessingLPT<nDim>>(postprocessingId_, data) {
  m_ppSolver = ppSolver_;
}

/**
 * \brief init arrays
 * \author Tim Wegmann
 * \date Jan. 2021
 *
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::initSprayData() {
  PostProcessing<nDim, PostProcessingLPT<nDim>>::initSprayData();
  TRACE();

  if(solver().domainId() == 0) {
    cerr << "Allocating Post-processing LPT data!" << endl;
  }

  // gatherd in the postProcTimeInterval
  mAlloc(m_particleCV, m_sprayDataSize, 2 + nDim, "m_particleCV", F0, AT_);
  mAlloc(m_particlePen, m_sprayDataSize, 19, "m_particlePen", F0, AT_);
  mAlloc(m_sprayStat, m_sprayDataSize, 5, "m_sprayStat", F0, AT_);
  mAlloc(m_sprayTimes, m_sprayDataSize, 3, "m_sprayTimes", F0, AT_);

  if(solver().m_sprayModel == nullptr) return;

  // gatherd every timeStep
  mAlloc(m_injectionData, this->m_sprayWriteInterval, 15, "m_injectionData", F0, AT_);
  mAlloc(m_cInjData, this->m_sprayWriteInterval + 1, nDim + 5, "m_cInjData", F0, AT_);

  if(!solver().isActive()) return;

  array<MFloat, nDim + 5> loadedVariables = {};
  // load the sum of the injected mass from last file
  if(solver().m_restartFile && solver().domainId() == 0) {
    ifstream readFile;
    stringstream filename;
    filename.clear();
    filename.str("");
    filename << solver().restartDir() << "injSums";
    readFile.open(filename.str(), ios_base::in);
    if(!readFile) {
      filename.clear();
      filename.str("");
      filename << solver().restartDir() << "injSums_" << globalTimeStep;
      readFile.open(filename.str(), ios_base::in);
      if(!readFile) {
        mTerm(1, AT_, "Error reading injection file!");
      }
    } else {
      // file found -> create backup
      MString backupName = solver().restartDir() + "injSums_" + to_string(globalTimeStep);
      std::ifstream src(filename.str(), std::ios::binary);
      std::ofstream dst(backupName, std::ios::binary);
      dst << src.rdbuf();
    }
    readFile.seekg(-1, ios_base::end);

    MChar ch;
    MString lastLine;
    readFile.get(ch);
    if(ch == '\n') {
      readFile.seekg(-2, ios::cur);
      readFile.seekg(-1, ios::cur);
      readFile.get(ch);
      while(ch != '\n') {
        readFile.seekg(-2, std::ios::cur);
        readFile.get(ch);
      }
      getline(readFile, lastLine);
    } else {
      mTerm(1, AT_, "Wrong injection file-format!");
    }
    stringstream ss(lastLine);
    MFloat timeStep;
    ss >> timeStep;

    if((MInt)timeStep == globalTimeStep) {
      // matching time-step -> load data
      for(MInt i = 0; i < nDim + 5; i++) {
        ss >> loadedVariables[i];
        cerr << "Loaded injected Data " << i << " is " << loadedVariables[i] << endl;
      }
    } else if((MInt)timeStep > globalTimeStep) {
      cerr << "Searching correct line in injection data-file for " << globalTimeStep << " from last output " << timeStep
           << endl;
      while(timeStep > globalTimeStep) {
        ss.str("");
        readFile.seekg(-2, std::ios::cur);
        readFile.seekg(-1, std::ios::cur);
        readFile.get(ch);
        const MInt oldPos = readFile.tellg();
        while(ch != '\n') {
          readFile.seekg(-2, std::ios::cur);
          readFile.get(ch);
        }
        MString lastLine2;
        getline(readFile, lastLine2);
        const MInt newPos = readFile.tellg();
        const MInt difPos = -(newPos - oldPos);
        readFile.seekg(difPos, std::ios::cur);
        ss.str(lastLine2);
        ss >> timeStep;
      }
      for(MInt i = 0; i < nDim + 5; i++) {
        ss >> loadedVariables[i];
        cerr << "Loaded injected Data " << i << " is " << loadedVariables[i] << endl;
      }

    } else if((MInt)timeStep < globalTimeStep) {
      for(MInt i = 0; i < nDim + 5; i++) {
        loadedVariables[i] = 0;
      }
      readFile.close();

      // if the restart timeStep is smaller than the last computed/written timeStep info
      //-> search for matching timeStep in the backup file
      filename.clear();
      filename.str("");
      filename << solver().restartDir() << "injSums_" << globalTimeStep;
      readFile.open(filename.str(), ios_base::in);
      if(readFile) {
        mTerm(1, AT_, "Error reading injection file!");
      }
      readFile.seekg(-1, ios_base::end);
      readFile.get(ch);
      MString lastLine2;
      if(ch == '\n') {
        readFile.seekg(-2, ios::cur);
        readFile.seekg(-1, ios::cur);
        readFile.get(ch);
        while(ch != '\n') {
          readFile.seekg(-2, std::ios::cur);
          readFile.get(ch);
        }
        getline(readFile, lastLine2);
      } else {
        mTerm(1, AT_, "Wrong injection file-format!");
      }
      stringstream ss2(lastLine2);
      ss2 >> timeStep;
      if((MInt)timeStep == globalTimeStep) {
        for(MInt i = 0; i < nDim + 5; i++) {
          ss2 >> loadedVariables[i];
          if(solver().domainId() == 0) {
            cerr << "Loaded injected Data " << i << " is " << loadedVariables[i] << endl;
          }
        }
      } else {
        cerr << "Correct Timestep count not be found at the end of in any of the injSums files!"
             << "Either trunkate the files of implement version which checks through all lines!" << endl;
      }
    }
    readFile.close();

    for(MInt i = 0; i < nDim + 5; i++) {
      m_cInjData[this->m_sprayWriteInterval][i] = loadedVariables[i];
    }
  }

  if(solver().m_restartFile) {
    MPI_Allreduce(MPI_IN_PLACE, &(m_cInjData[0][0]), (nDim + 5) * (this->m_sprayWriteInterval + 1), MPI_DOUBLE, MPI_MAX,
                  solver().mpiComm(), AT_, "INPLACE", "m_cInjData");
  }
}

/**
 * \brief computes conervation variables (mass/momentum/energy!)
 * \author Sven Berger, Tim Wegmann
 * \date Jan. 2021
 *
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::particleMass() {
  TRACE();

  // calculate conservative particle values
  MFloat partMass = 0;
  MFloat partE = 0;
  MFloat partMomentum[nDim] = {};
  for(MInt i = 0; i < nDim; i++) {
    partMomentum[i] = 0;
  }
  for(MInt id = 0; id < solver().a_noParticles(); ++id) {
    if(solver().m_partList[id].isInvalid()) continue;
    MFloat mass = solver().m_partList[id].sphericalMass() * solver().m_partList[id].m_noParticles;
    partMass += mass;
    partE +=
        solver().m_material->cp() / solver().m_material->gammaMinusOne() * mass * solver().m_partList[id].m_temperature;
    MFloat velMagSquared = 0;
    for(MInt i = 0; i < nDim; i++) {
      partMomentum[i] += (mass * solver().m_partList[id].m_velocity[i]);
      velMagSquared += POW2(solver().m_partList[id].m_velocity[i]);
    }
    partE += 0.5 * mass * velMagSquared;
  }
  m_particleCV[m_sprayDataStep][0] = partMass;
  m_particleCV[m_sprayDataStep][1] = partE;
  for(MInt i = 0; i < nDim; i++) {
    m_particleCV[m_sprayDataStep][2 + i] = partMomentum[i];
  }
}


/**
 * \brief calculate vertical and horizontal penetration from given coordinate
 * \author Sven Berger, Tim Wegmann
 * \date Jan. 2021
 *
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::particlePenetration() {
  TRACE();

  // different penetration measurements
  // 0: overal all maximum
  // 1: parcel with 90% liquid mass fraction
  // 2: parcel with LVFlimit
  // 3: parcel with LVFlimit on refLvl
  const MInt noPen = 4;
  array<array<MFloat, noPen>, nDim + 1> maxLiqPen{};

  const MFloat massFracLimit = 0.90;
  const MFloat LVFlimit = 0.025; // 0.01 as first estimate
  const MInt refLvl = 10;

  // spray-width measurements in certain plans
  static constexpr MInt nPlanes = 3;
  // static constexpr MFloat yPlaneBorder = 0.1/75.0;
  static constexpr MFloat yPlaneBorder = 0.001333;
  // static constexpr MFloat yPlanes[nPlanes] = {15.0/75.0 , 25.0/75.0, 35.0/75.0};
  static constexpr MFloat yPlanes[nPlanes] = {0.2, 0.3333, 0.46667};

  vector<MFloat> width(nPlanes);
  for(MInt n = 0; n < nPlanes; n++) {
    width[n] = 0.0;
  }

  auto assignLarger = [&](MFloat& A, MFloat b) {
    if(b > A) {
      A = b;
    }
  };


  for(MInt id = 0; id < solver().a_noParticles(); ++id) {
    array<MFloat, nDim + 1> distance{};
    for(MInt i = 0; i < nDim; i++) {
      distance[i] = abs(solver().m_partList[id].m_position[i] - solver().m_spawnCoord[i]);
    }

    distance[nDim] = sqrt(POW2(distance[0]) + POW2(distance[1]) + POW2(distance[2]));


    for(MInt i = 0; i < nDim + 1; i++) {
      assignLarger(maxLiqPen[i][0], distance[i]);
    }

    for(MInt i = 0; i < nPlanes; i++) {
      if(fabs(distance[1] - yPlanes[i]) < yPlaneBorder) {
        assignLarger(width[i], distance[0]);
        assignLarger(width[i], distance[2]);
      }
    }
    const MInt cellId = solver().m_partList[id].m_cellId;
    if(cellId < 0) continue;
    if(solver().a_volumeFraction(cellId) > LVFlimit) {
      for(MInt i = 0; i < nDim + 1; i++) {
        assignLarger(maxLiqPen[i][1], distance[i]);
      }
    }
  }

  for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
    const MFloat massFraction =
        solver().a_volumeFraction(cellId) * solver().m_material->density() / solver().a_fluidDensity(cellId);
    array<MFloat, nDim + 1> distance{};

    if(massFraction > massFracLimit) {
      array<MFloat, nDim> cellCoordinate{};
      for(MInt i = 0; i < nDim; i++) {
        cellCoordinate[i] = solver().c_coordinate(cellId, i);
        distance[i] = abs(cellCoordinate[i] - solver().m_spawnCoord[i]);
      }
      distance[nDim] = sqrt(POW2(distance[0]) + POW2(distance[1]) + POW2(distance[2]));

      for(MInt i = 0; i < nDim + 1; i++) {
        assignLarger(maxLiqPen[i][2], distance[i]);
      }
    }
    if(solver().a_level(cellId) == refLvl) {
      solver().reduceData(cellId, &solver().a_volumeFraction(0), 1);
      if(solver().a_volumeFraction(cellId) > LVFlimit) {
        array<MFloat, nDim> cellCoordinate{};
        for(MInt i = 0; i < nDim; i++) {
          cellCoordinate[i] = solver().c_coordinate(cellId, i);
          distance[i] = abs(cellCoordinate[i] - solver().m_spawnCoord[i]);
        }
        distance[nDim] = sqrt(POW2(distance[0]) + POW2(distance[1]) + POW2(distance[2]));

        for(MInt i = 0; i < nDim + 1; i++) {
          assignLarger(maxLiqPen[i][3], distance[i]);
        }
      }
    }
  }

  // save data:
  MInt it = 0;
  for(MInt n = 0; n < noPen; n++) {
    for(MInt i = 0; i < nDim + 1; i++) {
      m_particlePen[m_sprayDataStep][it++] = maxLiqPen[i][n];
    }
  }

  for(MInt n = 0; n < nPlanes; n++) {
    m_particlePen[m_sprayDataStep][it++] = width[n];
  }

  if(it > 19) {
    mTerm(1, AT_, "PP: Penetration-Storage not matching!");
  }
}

/**
 * \brief computes parcel statistics
 * \author Sven Berger, Tim Wegmann
 * \date Jan. 2021
 *
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::parcelStatistics() {
  TRACE();

  const MInt noPart = solver().a_noParticles();

  MLong noDroplets = 0;
  MFloat diameter = 0;
  MFloat surfaceArea = 0;
  MFloat volume = 0;
  for(MInt id = 0; id < noPart; ++id) {
    noDroplets += solver().m_partList[id].m_noParticles;
    diameter += solver().m_partList[id].m_noParticles * solver().m_partList[id].m_diameter;
    surfaceArea += solver().m_partList[id].m_noParticles * POW2(solver().m_partList[id].m_diameter);
    volume += solver().m_partList[id].m_noParticles * POW3(solver().m_partList[id].m_diameter);
  }


  // const MFloat SMD = volume / surfaceArea;
  // const MFloat meanD = diameter / noDroplets;

  m_sprayStat[m_sprayDataStep][0] = noPart;
  m_sprayStat[m_sprayDataStep][1] = noDroplets;
  m_sprayStat[m_sprayDataStep][2] = diameter;
  m_sprayStat[m_sprayDataStep][3] = surfaceArea;
  m_sprayStat[m_sprayDataStep][4] = volume;
}


template <MInt nDim>
MInt PostProcessingLPT<nDim>::getInjectionData() {
  TRACE();


  if(solver().m_sprayModel == nullptr) {
    return 0;
  }

  if(!solver().isActive()) {
    solver().m_injData.clear();
    return 0;
  }

  const MInt injectorCellId = solver().injectorCellId();
  MInt count = 0;

  // gather and sum all spray injection solver data
  for(auto it = solver().m_injData.begin(); it != solver().m_injData.end(); it++) {
    const MInt timeStep = it->first;
    if(solver().domainId() == 0) {
      m_injectionData[count][0] = timeStep;
    } else {
      m_injectionData[count][0] = 0.0;
    }

    for(MInt i = 0; i < 14; i++) {
      const MFloat data = (it->second)[i];
      m_injectionData[count][i + 1] = data;
    }
    count++;
  }
  count = 0;

  MPI_Allreduce(MPI_IN_PLACE, &(m_injectionData[0][0]), 15 * this->m_sprayWriteInterval, MPI_DOUBLE, MPI_SUM,
                solver().mpiComm(), AT_, "INPLACE", "m_injectionData");

  MPI_Allreduce(MPI_IN_PLACE, &(m_cInjData[0][0]), (nDim + 5) * (this->m_sprayWriteInterval + 1), MPI_DOUBLE, MPI_MAX,
                solver().mpiComm(), AT_, "INPLACE", "m_cInjData");

  if(injectorCellId > -1) {
    // get storage value
    for(MInt i = 0; i < nDim + 5; i++) {
      m_cInjData[0][i] = m_cInjData[this->m_sprayWriteInterval][i];
    }

    for(auto it = solver().m_injData.begin(); it != solver().m_injData.end(); it++) {
      // add data from:
      // injMass, injXMom, injYMom, injZMom, injEnergy, sumEvapMass, noRT-breakup, noKH-breakup
      for(MInt i = 0; i < nDim + 5; i++) {
        m_cInjData[count][i] += m_injectionData[count][7 + i];
      }
      count++;

      // set pervious mass as mass in the next timeStep
      for(MInt i = 0; i < nDim + 5; i++) {
        m_cInjData[count][i] = m_cInjData[count - 1][i];
      }
    }
    // put last value into storage
    for(MInt i = 0; i < nDim + 5; i++) {
      m_cInjData[this->m_sprayWriteInterval][i] = m_cInjData[count][i];
    }

  } else {
    for(MInt id = 0; id < this->m_sprayWriteInterval + 1; id++) {
      for(MInt j = 0; j < nDim + 5; j++) {
        m_cInjData[id][j] = 0;
      }
    }
  }

  solver().m_injData.clear();

  return count;
}

/** \brief init function for LPT particle solution file
 * \author Tim Wegmann
 * \date Jan 2021
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::initLPTSolutionFile() {
  TRACE();
  m_LPTSolutionInterval =
      Context::getSolverProperty<MInt>("particleSolutionInterval", m_postprocessingId, AT_, &m_LPTSolutionInterval);

  m_writeParLog = Context::getSolverProperty<MBool>("writeParLog", m_postprocessingId, AT_, &m_writeParLog);
  if(m_writeParLog) initParticleLog();
}

/**
 * \brief base function which is called every time-Step in postprocessInSolve
 * \author Tim Wegmann, Johannes Grafen
 * \date Jan. 2021
 *
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::writeLPTSolutionFile() {
  TRACE();

  if(!solver().isActive()) return;

  solver().crankAngleSolutionOutput();

  if(m_writeParLog && globalTimeStep % m_parLogInterval == 0) writeParticleLog();

  if((m_LPTSolutionInterval > 0 && globalTimeStep % m_LPTSolutionInterval == 0) || this->m_finalTimeStep
     || this->m_forceOutput) {
    if(solver().domainId() == 0) {
      cerr << "Writing particle File @" << globalTimeStep << "(" << solver().m_time << ")... ";
    }
    solver().writePartData();

    if(solver().m_collisions > 0) {
      solver().writeCollData();
    }

    if(solver().domainId() == 0) {
      cerr << " done. " << endl;
    }
  }
}

/**
 * \brief initialize particle.log
 * \author Johannes Grafen
 */
template <MInt nDim>
void PostProcessingLPT<nDim>::initParticleLog() {
  if(solver().m_restart && m_writeParLog) m_parLogApp = true;

  m_parLogInterval =
      Context::getSolverProperty<MInt>("particleLogInterval", m_postprocessingId, AT_, &m_parLogInterval);

  if(solver().domainId() == 0) {
    std::ofstream parLog;
    if(m_parLogApp) {
      parLog.open("particle.log", ios::app);
    } else {
      parLog.open("particle.log");
      const MString columns[9] = {"x", "y", "z", "vel_x", "vel_y", "vel_z", "Re", "C_D"};
      parLog << std::setw(5) << "t"
             << "\t";
      for(auto& i : columns) {
        parLog << std::setw(15) << i << "\t";
      }
    }
    parLog << std::endl;
    parLog.close();
  }
}

/**
 * \brief  write particle Log file, containing particle position and velocity
 * at timestep t
 *
 * \author Johannes Grafen
 */
template <MInt nDim>
void PostProcessingLPT<nDim>::writeParticleLog() {
  if(solver().globalNoParticles() > 1) {
    m_writeParLog = false;
    std::cerr << "writing particle.log is currently only implemented for single particles!" << std::endl;
    return;
  }

  if(solver().m_partList.size() < 1) return;


  std::ofstream parLog;
  LPTSpherical<nDim>& part = solver().m_partList[0];
  if(m_writeParLog) {
    parLog.open("particle.log", std::ios::app);

    parLog << std::setw(5) << globalTimeStep << "\t";
    for(const MFloat& i : part.m_position) {
      parLog << std::setw(10) << std::setprecision(12) << i << "\t";
    }
    for(const MFloat& i : part.m_velocity) {
      parLog << std::setw(10) << i << "\t";
    }

    const MFloat T = solver().a_fluidTemperature(part.m_cellId);
    const MFloat fluidViscosity = solver().m_material->dynViscosityFun(T);
    const MFloat relVel = part.magRelVel(&part.m_oldFluidVel[0], &part.m_oldVel[0]);

    const MFloat ReP = part.particleRe(relVel, part.m_oldFluidDensity, fluidViscosity) * part.s_Re;
    const MFloat CD = part.dragFactor(ReP) * 24 / ReP;

    parLog << std::setw(10) << ReP << "\t";
    parLog << std::setw(10) << CD << "\t" << std::endl;
    parLog.close();
  }
}

/** \brief init function for LPT particle statistics file
 * \author Johannes Grafen
 * \date Jul 2022
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::initParticleStatistics() {
  TRACE();

  if(solver().domainId() == 0) {
    std::ofstream Plog;
    Plog.open("ParticleStatistics.log", ios::app);
    Plog << std::setw(5) << "t"
         << "\t" << std::setw(10) << "VpRMS"
         << "\t" << std::setw(10) << "AverageReP"
         << "\t" << std::setw(10) << "EkinP" << endl;
    Plog.close();
  }
}
/**
 * \brief compute average quantites of particle phase
 *    average particle Reynolds number ReP
 *    the root mean square particle velocity vrms
 *    the total kinetic energy of the particles
 * \author Johannes Grafen
 * \date Jul 2022
 *
 **/
template <MInt nDim>
void PostProcessingLPT<nDim>::computeParticleStatistics() {
  TRACE();

  if((m_LPTSolutionInterval > 0 && globalTimeStep % m_LPTSolutionInterval == 0) || this->m_finalTimeStep
     || this->m_forceOutput) {
    const MInt globalNoPart = solver().globalNoParticles();
    MFloat RePSum = F0;
    MFloat VelmagnitudeSqSum = F0;
    MFloat EkinSum = F0;

    for(LPTSpherical<nDim>& part : solver().m_partList) {
      const MFloat T = solver().a_fluidTemperature(part.m_cellId);
      const MFloat fluidViscosity = solver().m_material->dynViscosityFun(T);
      const MFloat relVel = part.magRelVel(&part.m_oldFluidVel[0], &part.m_oldVel[0]);

      RePSum += part.particleRe(relVel, part.m_oldFluidDensity, fluidViscosity) * part.s_Re;

      MFloat VelmagnitudeSq = F0;
      for(MInt j = 0; j < nDim; j++) {
        VelmagnitudeSq += POW2(part.m_velocity[j]);
      }
      VelmagnitudeSqSum += VelmagnitudeSq;

      EkinSum += part.sphericalMass() * VelmagnitudeSq;
    }

    MPI_Allreduce(MPI_IN_PLACE, &RePSum, 1, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_, "INPLACE", "RePsum");

    MFloat AverageReP = RePSum / globalNoPart;

    MPI_Allreduce(MPI_IN_PLACE, &VelmagnitudeSqSum, 1, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_, "INPLACE",
                  "VelMagnitudesum");

    MFloat VpRMS = sqrt(VelmagnitudeSqSum / globalNoPart);

    MPI_Allreduce(MPI_IN_PLACE, &EkinSum, 1, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_, "INPLACE", "EkinSum");

    MFloat EkinP = 0.5 * EkinSum;

    if(solver().domainId() == 0) {
      std::ofstream Plog;
      Plog.open("ParticleStatistics.log", ios::app);
      Plog << std::setw(5) << globalTimeStep << "\t" << std::setw(10) << VpRMS << "\t" << std::setw(10) << AverageReP
           << "\t" << std::setw(10) << EkinP << endl;
      Plog.close();
    }
  } else {
    return;
  }
}


template class PostProcessingLPT<3>;
