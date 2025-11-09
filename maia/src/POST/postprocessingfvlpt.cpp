// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessingfvlpt.h"
#include "IO/cpptoml.h"
#include "postdata.h"
#include "postprocessing.cpp"


using namespace std;

template <MInt nDim, class SysEqn>
PostProcessingFvLPT<nDim, SysEqn>::PostProcessingFvLPT(MInt postprocessingId_,
                                                       PostData<nDim>* data,
                                                       FvCartesianSolverXD<nDim, SysEqn>* ppFvSolver_,
                                                       LPT<nDim>* ppLptSolver_)
  : PostProcessingInterface(postprocessingId_),
    PostProcessingLPT<nDim>(postprocessingId_, data, ppLptSolver_),
    PostProcessingFv<nDim, SysEqn>(postprocessingId_, data, ppFvSolver_) {
  m_ppSolverLpt = ppLptSolver_;
  m_ppSolverFv = ppFvSolver_;

  m_fvLPTSpeciesId = 0;
  m_fvLPTSpeciesId = Context::getSolverProperty<MInt>("fvLptSpeciesId", lpt().solverId(), AT_, &m_fvLPTSpeciesId);
}

/**
 * \brief base function which is called every time-Step in postprocessInSolve
 * \author Tim Wegmann
 * \date Jan. 2021
 *
 **/
template <MInt nDim, class SysEqn>
void PostProcessingFvLPT<nDim, SysEqn>::computeSprayData() {
  TRACE();


  if(globalTimeStep % this->m_sprayComputeInterval == 0) {
    updateData();
    advanceDataStep();
  }


  if(globalTimeStep % this->m_sprayWriteInterval == 0) {
    writeSprayData();
    resetDataStep();
  }
}

/**
 * \brief updates all post-processing data from postprocessing-fv or postprocessing-lpt
 * \author Tim Wegmann
 * \date Jan. 2021
 *
 **/
template <MInt nDim, class SysEqn>
void PostProcessingFvLPT<nDim, SysEqn>::updateData() {
  TRACE();


  m_sprayTimes[m_sprayDataStep][0] = globalTimeStep;
  m_sprayTimes[m_sprayDataStep][1] = fvSolver().physicalTime();
  m_sprayTimes[m_sprayDataStep][2] = fvSolver().crankAngle(fvSolver().physicalTime(), 0);

  fvSolver().stopLoadTimer(AT_);
  if(lpt().isActive()) {
    lpt().startLoadTimer(AT_);
    this->particleMass();
    this->particlePenetration();
    this->parcelStatistics();
    lpt().stopLoadTimer(AT_);
  }

  fvSolver().startLoadTimer(AT_);

  if(fvSolver().isActive()) {
    this->vapourMass(m_fvLPTSpeciesId);
    this->vapourPenetration(lpt().m_spawnCoord);
  }
}

template <MInt nDim, class SysEqn>
void PostProcessingFvLPT<nDim, SysEqn>::writeSprayData() {
  TRACE();

  static MBool firstTime = true;

  if(this->m_finalTimeStep && fvSolver().domainId() == 0) {
    cerr << "Writing final postProcessing Data" << endl;
  }

  if(firstTime && fvSolver().domainId() == 0) {
    cerr0 << "Writing first postProcessing Data" << endl;
  }

  fvSolver().stopLoadTimer(AT_);

  const MBool timerEnabled = fvSolver().dlbTimersEnabled();

  if(timerEnabled) {
    maia::dlb::g_dlbTimerController.disableAllDlbTimers();
  }

  const MBool hasInjection = lpt().m_sprayModel != nullptr ? true : false;
  MInt injDataSize = this->getInjectionData();
  MPI_Allreduce(MPI_IN_PLACE, &injDataSize, 1, MPI_INT, MPI_MAX, fvSolver().mpiComm(), AT_, "INPLACE", "injDataSize");
  const MInt dataSize = (this->m_finalTimeStep || firstTime) ? m_sprayDataStep : m_sprayDataSize;
  const MInt vapourPenSize = nDim == 3 ? 16 : 12;

  // rank with injector writes data
  const MInt injectorCellId = lpt().isActive() ? lpt().injectorCellId() : -1;

  if(hasInjection && injDataSize != dataSize) {
    cerr0 << "Spray-PP data differs: " << injDataSize << " inj-Size " << dataSize << " sprayData-Size" << endl;
  }

  if(dataSize > 0) {
    // exchange data:
    MPI_Allreduce(MPI_IN_PLACE, &m_particleCV[0][0], m_sprayDataSize * (2 + nDim), MPI_DOUBLE, MPI_SUM,
                  fvSolver().mpiComm(), AT_, "INPLACE", "m_particleCV");

    MPI_Allreduce(MPI_IN_PLACE, &m_sprayStat[0][0], m_sprayDataSize * 5, MPI_DOUBLE, MPI_SUM, fvSolver().mpiComm(), AT_,
                  "INPLACE", "m_sprayStat");

    MPI_Allreduce(MPI_IN_PLACE, &m_particlePen[0][0], m_sprayDataSize * 19, MPI_DOUBLE, MPI_MAX, fvSolver().mpiComm(),
                  AT_, "INPLACE", "m_particlePen");

    MPI_Allreduce(MPI_IN_PLACE, &m_vapourCV[0][0], m_sprayDataSize * (5 + 2 * nDim), MPI_DOUBLE, MPI_SUM,
                  fvSolver().mpiComm(), AT_, "INPLACE", "m_vapourCV");

    MPI_Allreduce(MPI_IN_PLACE, &m_vapourPen[0][0], m_sprayDataSize * vapourPenSize, MPI_DOUBLE, MPI_MAX,
                  fvSolver().mpiComm(), AT_, "INPLACE", "m_vapourPen");
  }

  if(injectorCellId >= 0) {

    ofstream ofl;

    // backup a few older files
    if(firstTime) {
      // file with conservation variables
      MString fileNameCheck = lpt().outputDir() + "sprayCons";
      if(fileExists(fileNameCheck)) {
        const MString fileName = lpt().outputDir() + "sprayCons_" + to_string(lpt().m_restartTimeStep);
        if(!fileExists(fileName)) copyFile(fileNameCheck, fileName);
      }
      if(!lpt().m_restartFile) {
        ofl.open(fileNameCheck, ios_base::out | ios_base::trunc);
        ofl << "#0:ts 1:time 2:CAD 3:pMass 4:pEnergy 5-7:pMomentum 8:vMass 9:aMass 10:gMass 11:gEnergy 12:gMomX "
               "13:gMomY 14:gMomZ";
        ofl << endl;
        ofl.close();
      }

      // file with liquid penetration lengths
      fileNameCheck = lpt().outputDir() + "liquidPen";
      if(fileExists(fileNameCheck)) {
        const MString fileName = lpt().outputDir() + "liquidPen_" + to_string(lpt().m_restartTimeStep);
        if(!fileExists(fileName)) copyFile(fileNameCheck, fileName);
      }
      if(!lpt().m_restartFile) {
        ofl.open(fileNameCheck, ios_base::out | ios_base::trunc);
        ofl << "#0:ts 1:time 2:CAD 3-5:liquidPen"
               "6-8:liquidWidth";
        ofl << endl;
        ofl.close();
      }

      // file with vapour penetration lengths
      fileNameCheck = lpt().outputDir() + "vapourPen";
      if(fileExists(fileNameCheck)) {
        const MString fileName = lpt().outputDir() + "vapourPen_" + to_string(lpt().m_restartTimeStep);
        if(!fileExists(fileName)) copyFile(fileNameCheck, fileName);
      }
      if(!lpt().m_restartFile) {
        ofl.open(fileNameCheck, ios_base::out | ios_base::trunc);
        ofl << "#0:ts 1:time 2:CAD 3-5:v01Pen 6-8:v1Pen 9-11:v5Pen";
        ofl << endl;
        ofl.close();
      }

      // file with spray statistics
      fileNameCheck = lpt().outputDir() + "sprayStats";
      if(fileExists(fileNameCheck)) {
        const MString fileName = lpt().outputDir() + "sprayStats_" + to_string(lpt().m_restartTimeStep);
        if(!fileExists(fileName)) copyFile(fileNameCheck, fileName);
      }

      if(!lpt().m_restartFile) {
        ofl.open(fileNameCheck, ios_base::out | ios_base::trunc);
        ofl << "#0:ts 1:time 2:CAD 3:noParcles 4:noDroplets 5:SMS 6:avgDiameter";
        ofl << endl;
        ofl.close();
      }

      // file with injection statistics
      if(hasInjection) {
        fileNameCheck = lpt().outputDir() + "injStats";
        if(fileExists(fileNameCheck)) {
          const MString fileName = lpt().outputDir() + "injStats_" + to_string(lpt().m_restartTimeStep);
          if(!fileExists(fileName)) copyFile(fileNameCheck, fileName);
        }
        if(!lpt().m_restartFile) {
          ofl.open(fileNameCheck, ios_base::out | ios_base::trunc);
          ofl << "#0:ts 1:timeSinceSOI 2:injProgress 3:rampFactor 4:noInjPart 5:noInjDrop 6:injDiameter 7:injMass "
                 "8:injXMom 9:injYMom 10:injZMom 11:injEnergy";
          ofl << endl;
          ofl.close();
        }

        // file with sum of conservative variables from injection
        fileNameCheck = lpt().outputDir() + "injSums";
        if(fileExists(fileNameCheck)) {
          const MString fileName = lpt().outputDir() + "injSums_" + to_string(lpt().m_restartTimeStep);
          if(!fileExists(fileName)) copyFile(fileNameCheck, fileName);
        }
        if(!lpt().m_restartFile) {
          ofl.open(fileNameCheck, ios_base::out | ios_base::trunc);
          ofl << "#0:ts 1:injMass 2-4:injMomentum 5:injEnergy 6:evapMass";
          ofl << endl;
          ofl.close();
        }
      }
    }


    MString fileName = lpt().outputDir() + "sprayCons";
    ofl.open(fileName, ios_base::out | ios_base::app);

    if(ofl.is_open() && ofl.good()) {
      for(MInt i = 0; i < dataSize; i++) {
        ofl << m_sprayTimes[i][0] << setprecision(7) << " " << m_sprayTimes[i][1] << " " << m_sprayTimes[i][2] << " "
            << setprecision(12);
        for(MInt v = 0; v < 2 + nDim; v++) {
          ofl << m_particleCV[i][v] << " ";
        }
        for(MInt v = 0; v < 5 + 2 * nDim; v++) {
          ofl << m_vapourCV[i][v] << " ";
        }
        ofl << endl;
      }
    }
    ofl.close();

    fileName = lpt().outputDir() + "liquidPen";
    ofl.open(fileName, ios_base::out | ios_base::app);
    if(ofl.is_open() && ofl.good()) {
      for(MInt i = 0; i < dataSize; i++) {
        ofl << m_sprayTimes[i][0] << setprecision(7) << " " << m_sprayTimes[i][1] << " " << setprecision(7)
            << m_sprayTimes[i][2] << " " << setprecision(8);
        for(MInt j = 0; j < 19; j++) {
          ofl << m_particlePen[i][j] << " ";
        }
        ofl << endl;
      }
    }
    ofl.close();

    fileName = lpt().outputDir() + "vapourPen";
    ofl.open(fileName, ios_base::out | ios_base::app);
    if(ofl.is_open() && ofl.good()) {
      for(MInt i = 0; i < dataSize; i++) {
        ofl << m_sprayTimes[i][0] << setprecision(7) << " " << m_sprayTimes[i][1] << " " << setprecision(7)
            << m_sprayTimes[i][2] << " " << setprecision(8);
        for(MInt j = 0; j < vapourPenSize; j++) {
          ofl << m_vapourPen[i][j] << " ";
        }
        ofl << endl;
      }
    }
    ofl.close();

    fileName = lpt().outputDir() + "sprayStats";
    ofl.open(fileName, ios_base::out | ios_base::app);
    if(ofl.is_open() && ofl.good()) {
      for(MInt i = 0; i < dataSize; i++) {
        ofl << m_sprayTimes[i][0] << setprecision(7) << " " << m_sprayTimes[i][1] << " " << m_sprayTimes[i][2] << " "
            << setprecision(12);
        for(MInt j = 0; j < 5; j++) {
          ofl << m_sprayStat[i][j] << " ";
        }
        ofl << endl;
      }
    }
    ofl.close();

    if(hasInjection) {
      fileName = lpt().outputDir() + "injStats";
      ofl.open(fileName, ios_base::out | ios_base::app);
      if(ofl.is_open() && ofl.good()) {
        for(MInt i = 0; i < injDataSize; i++) {
          ofl << setprecision(12);
          for(MInt j = 0; j < 15; j++) {
            ofl << m_injectionData[i][j] << " ";
          }
          ofl << endl;
        }
      }
      ofl.close();

      fileName = lpt().outputDir() + "injSums";
      ofl.open(fileName, ios_base::out | ios_base::app);
      if(ofl.is_open() && ofl.good()) {
        for(MInt i = 0; i < injDataSize; i++) {
          ofl << m_injectionData[i][0] << " " << setprecision(12);
          for(MInt j = 0; j < nDim + 5; j++) {
            ofl << m_cInjData[i][j] << " ";
          }
          ofl << endl;
        }
      }
      ofl.close();
    }
  }
  firstTime = false;

  // reset data
  for(MInt i = 0; i < m_sprayDataSize; i++) {
    for(MInt j = 0; j < 5; j++) {
      m_sprayStat[i][j] = 0;
    }
    for(MInt j = 0; j < 2 + nDim; j++) {
      m_particleCV[i][j] = 0;
    }
    for(MInt j = 0; j < 19; j++) {
      m_particlePen[i][j] = 0;
    }
    if(hasInjection) {
      for(MInt j = 0; j < 15; j++) {
        m_injectionData[i][j] = 0;
      }
    }
  }

  if(timerEnabled) {
    maia::dlb::g_dlbTimerController.enableAllDlbTimers();
  }
  fvSolver().startLoadTimer(AT_);
}

template class PostProcessingFvLPT<3, FvSysEqnNS<3>>;
template class PostProcessingFvLPT<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class PostProcessingFvLPT<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class PostProcessingFvLPT<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
