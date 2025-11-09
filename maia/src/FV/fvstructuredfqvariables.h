// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDFQVARIABLES_H
#define STRUCTUREDFQVARIABLES_H

class StructuredFQVariables {
 public:
  StructuredFQVariables() {
    // this is the start configuration,
    // depending on which fq fields are needed
    // it will be altered later
    AVG_U = 0;
    AVG_V = 1;
    AVG_W = 2;
    AVG_P = 3;
    AVG_RHO = 4;
    NU_T = 5;
    SPONGE_RHO = 6;
    SPONGE_RHO_E = 7;
    PRESSURE = 8;
    LAMBDA2 = 9;
    VORTX = 10;
    VORTY = 11;
    VORTZ = 12;
    VORTICITY = new MInt[3];
    VORTICITY[0] = VORTX;
    VORTICITY[1] = VORTY;
    VORTICITY[2] = VORTZ;
    CELLID = 13;
    BLOCKID = 14;
    SPONGE_FACTOR = 15;
    SLOPEX = 16;
    SLOPEY = 17;
    SLOPEZ = 18;
    SLOPE = new MInt[3];
    SLOPE[0] = SLOPEX;
    SLOPE[1] = SLOPEY;
    SLOPE[2] = SLOPEZ;
    MU_T = 19;
    WALLDISTANCE = 20;
    MU_L = 21;
    UTAU = 22;
    // Porous
    POROSITY = 23;
    DARCY = 24;
    FORCH = 25;
    NORMALX = 26;
    NORMALY = 27;
    NORMALZ = 28;
    NORMAL = new MInt[3];
    NORMAL[0] = NORMALX;
    NORMAL[1] = NORMALY;
    NORMAL[2] = NORMALZ;
    LIMITERVISC = 29;
    POROUS_INDICATOR = 30;
    UTAU2 = 31;

    noFQVariables = 0;
    maxNoFQVariables = 32;
    noFQBoxOutput = 0;
  }

  ~StructuredFQVariables() {
    delete[] VORTICITY;
    delete[] SLOPE;
    delete[] NORMAL;
  }

  MInt AVG_U;
  MInt AVG_V;
  MInt AVG_W;
  MInt AVG_P;
  MInt AVG_RHO;
  MInt NU_T;
  MInt SPONGE_RHO;
  MInt SPONGE_RHO_E;
  MInt PRESSURE;
  MInt LAMBDA2;
  MInt VORTX;
  MInt VORTY;
  MInt VORTZ;
  MInt BLOCKID;
  MInt CELLID;
  MInt SPONGE_FACTOR;
  MInt* VORTICITY;
  MInt SLOPEX;
  MInt SLOPEY;
  MInt SLOPEZ;
  MInt* SLOPE;
  MInt MU_T;
  MInt WALLDISTANCE;
  MInt MU_L;
  MInt UTAU;
  MInt UTAU2;
  MInt POROSITY;
  MInt DARCY;
  MInt FORCH;
  MInt NORMALX;
  MInt NORMALY;
  MInt NORMALZ;
  MInt* NORMAL;
  MInt LIMITERVISC;
  MInt POROUS_INDICATOR;

  MInt noFQVariables;
  MInt noFQBoxOutput;
  MInt maxNoFQVariables;

  MInt* neededFQVariables = nullptr;
  MBool* outputFQVariables = nullptr; // write output only for certain FQ variables
  MBool* boxOutputFQVariables = nullptr;
  MBool* loadedFromRestartFile = nullptr;
  std::vector<MString> fqNames;
  std::vector<MBool> fqWriteOutput;
  std::vector<MBool> fqWriteOutputBoxes;

  void activateFQField(MInt fqIndex, MInt currentPos, MBool writeOutput, MBool writeOutputBoxes) {
    fqWriteOutput.push_back(writeOutput);
    fqWriteOutputBoxes.push_back(writeOutputBoxes);

    switch(fqIndex) {
      case 0:
        AVG_U = currentPos;
        fqNames.push_back("avgU");
        m_log << "FQ-Field: AVG_U activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 1:
        AVG_V = currentPos;
        fqNames.push_back("avgV");
        m_log << "FQ-Field: AVG_V activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 2:
        AVG_W = currentPos;
        fqNames.push_back("avgW");
        m_log << "FQ-Field: AVG_P activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 3:
        AVG_P = currentPos;
        fqNames.push_back("avgP");
        m_log << "FQ-Field: AVG_P activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 4:
        AVG_RHO = currentPos;
        fqNames.push_back("avgRho");
        m_log << "FQ-Field: AVG_RHO activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 5:
        NU_T = currentPos;
        fqNames.push_back("nu_t");
        m_log << "FQ-Field: NU_T activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 6:
        SPONGE_RHO = currentPos;
        fqNames.push_back("spongeRho");
        m_log << "FQ-Field: SPONGE_RHO activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 7:
        SPONGE_RHO_E = currentPos;
        fqNames.push_back("spongeRhoE");
        m_log << "FQ-Field: SPONGE_RHO_E activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 8:
        PRESSURE = currentPos;
        fqNames.push_back("pressure");
        m_log << "FQ-Field: PRESSURE activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 9:
        LAMBDA2 = currentPos;
        fqNames.push_back("lambda2");
        m_log << "FQ-Field: LAMBDA_2 activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 10:
        VORTX = currentPos;
        VORTICITY[0] = currentPos;
        fqNames.push_back("vortX");
        m_log << "FQ-Field: x-Vorticity activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 11:
        VORTY = currentPos;
        VORTICITY[1] = currentPos;
        fqNames.push_back("vortY");
        m_log << "FQ-Field: y-Vorticity activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 12:
        VORTZ = currentPos;
        VORTICITY[2] = currentPos;
        fqNames.push_back("vortZ");
        m_log << "FQ-Field: z-Vorticity activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 13:
        CELLID = currentPos;
        fqNames.push_back("cellId");
        m_log << "FQ-Field: cellId activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 14:
        BLOCKID = currentPos;
        fqNames.push_back("blockId");
        m_log << "FQ-Field: solverId activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 15:
        SPONGE_FACTOR = currentPos;
        fqNames.push_back("spongeFactor");
        m_log << "FQ-Field: SPONGE_FACTOR activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 16:
        SLOPEX = currentPos;
        SLOPE[0] = currentPos;
        fqNames.push_back("slopeX");
        m_log << "FQ-Field: x-Slope activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 17:
        SLOPEY = currentPos;
        SLOPE[1] = currentPos;
        fqNames.push_back("slopeY");
        m_log << "FQ-Field: y-Slope activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 18:
        SLOPEZ = currentPos;
        SLOPE[2] = currentPos;
        fqNames.push_back("slopeZ");
        m_log << "FQ-Field: z-Slope activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 19:
        MU_T = currentPos;
        fqNames.push_back("mu_t");
        m_log << "FQ-Field: mu_t activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 20:
        WALLDISTANCE = currentPos;
        fqNames.push_back("wallDistance");
        m_log << "FQ-Field: wallDistance activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 21:
        MU_L = currentPos;
        fqNames.push_back("mu_l");
        m_log << "FQ-Field: mu_l activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 22:
        UTAU = currentPos;
        fqNames.push_back("utau");
        m_log << "FQ-Field: UTAU activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 23:
        POROSITY = currentPos;
        fqNames.push_back("porosity");
        m_log << "FQ-Field: POROSITY activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 24:
        DARCY = currentPos;
        fqNames.push_back("darcy");
        m_log << "FQ-Field: DARCY activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 25:
        FORCH = currentPos;
        fqNames.push_back("forch");
        m_log << "FQ-Field: FORCH activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 26:
        NORMALX = currentPos;
        NORMAL[0] = currentPos;
        fqNames.push_back("normalx");
        m_log << "FQ-Field: x-NORMAL activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 27:
        NORMALY = currentPos;
        NORMAL[1] = currentPos;
        fqNames.push_back("normaly");
        m_log << "FQ-Field: y-NORMAL activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 28:
        NORMALZ = currentPos;
        NORMAL[2] = currentPos;
        fqNames.push_back("normalz");
        m_log << "FQ-Field: z-NORMAL activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 29:
        LIMITERVISC = currentPos;
        fqNames.push_back("limitervisc");
        m_log << "FQ-Field: LIMITERVISC activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 30:
        POROUS_INDICATOR = currentPos;
        fqNames.push_back("porous_indicator");
        m_log << "FQ-Field: POROUS_INDICATOR activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      case 31:
        UTAU2 = currentPos;
        fqNames.push_back("utau2");
        m_log << "FQ-Field: UTAU2 activated - Write Output: " << fqWriteOutput[fqWriteOutput.size() - 1]
              << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size() - 1] << std::endl;
        break;
      default:
        mTerm(1, AT_, "Can't find the FQ field!");
        break;
    }
  }
};

#endif
