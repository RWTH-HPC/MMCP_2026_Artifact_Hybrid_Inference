// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CARTESIANNETCDF_H
#define CARTESIANNETCDF_H

#include "INCLUDE/maiatypes.h"

/*! \class CartesianNetcdf
    \date begin: 00.00.00 change: 00.00.00

    \brief define the names of all variables and attributes in the netcdf file
    \details NOTE: The corresponding .cpp file (cartesiannetcdf.cpp) was removed in a
             commit to SVN since it was not compiled and not used.
 */
class CartesianNetcdf {
 public:
  CartesianNetcdf()
    : nDim("nDim"),
      noCells("noCells"),
      noConsistentCells("noConsistentCells"),
      domainId("domainId"),
      adjacentDomains("adjacentDomains"),
      windowCellCode("windowCellCode"),
      originalDomain("originalDomain"),
      currentDomain("currentDomain"),
      originalId("originalId"),
      periodicDomain("periodicDomain"),
      quasiPeriodicDomain("quasiPeriodicDomain"),
      noChildIds("noChildIds"),
      noCellsDirections("noCellsDirections"),
      noNghbrIds("noNghbrIds"),
      lengthLevel0("lengthLevel0"),
      maxLevel("maxLevel"),
      parentIds("parentIds"),
      noCellChildIds("noCellChildIds"),
      childIds("childIds"),
      noCellNghbrIds("noCellNghbrIds"),
      nghbrIds("nghbrIds"),
      intCoordinates("intCoordinates"),
      level("level"),
      per("per"),
      qper("qper"),

      noBndCnds("noBndCnds"),
      bndCndDim("bndCndDim."),
      bndCndIds("bndCndIds"),
      bndCellIds("bndCellIds"),
      bndCellSideIds("bndCellSideIds"),

      coordinates("coordinates"),
      bodyCoordinates("bodyCoordinates"),

      cellIds("cellIds"),
      volume("volume"),
      adaptation("adaptation"),
      cutSideIds("cutSideIds"),
      noCutSideIds("noCutSideIds"),
      nonFluidSideIds("nonFluidSideIds"),
      noNonFluidSideIds("noNonFluidSideIds"),
      bndryCndId("bndryCndId"),
      slope("slope"),
      curvature("curvature"),
      area("area"),
      normalVctr("normalVctr"),
      cutCoordinates("cutCoordinates"),
      noCutPoints("noCutPoints"),
      cutWedgeIds("cutWedgeIds"),
      segmentIds("segmentIds"),
      bodyIds("bodyIds"),
      srfcs("srfcs"),
      noSrfcs("noSrfcs"),


      rightHandSide("rightHandSide"),
      rcnstrctnConstants("rcnstrctnConstants"),
      noRcnstrctnNghbrIds("noRcnstrctnNghbrIds"),
      rcnstrctnNghbrIds("rcnstrctnNghbrIds"),
      u_Velocity("u_Velocity"),
      v_Velocity("v_Velocity"),
      w_Velocity("w_Velocity"),
      temperature("temperature"),
      energy("energy"),
      rho("rho"),
      passiveScalar("passiveScalar"),
      minValues("minValues"),
      maxValues("maxValues"),
      noSamples("noSamples"),
      globalTimeStep("globalTimeStep"),
      del("del"),
      time("time"),
      physicalTime("physicalTime"),
      firstMaxResidual("firstMaxResidual"),
      firstAvrgResidual("firstAvrgResidual"),
      levelSetFunction("levelSetFunction"),
      oldLevelSetFunction("oldLevelSetFunction"),
      extensionVelocity("extensionVelocity"),
      semiLagrangeXShiftRef("semiLagrangeXShiftRef"),

      cellIsInactive("cellIsInactive"),

      flowCellId("flowCellId"),
      bndryId("bndryId") {
    variables[0] = "variables0";
    variables[1] = "variables1";
    variables[2] = "variables2";
    variables[3] = "variables3";
    variables[4] = "variables4";
    variables[5] = "variables5";

    oldVariables[0] = "oldVariables0";
    oldVariables[1] = "oldVariables1";
    oldVariables[2] = "oldVariables2";
    oldVariables[3] = "oldVariables3";
    oldVariables[4] = "oldVariables4";
    oldVariables[5] = "oldVariables5";

    variables[6] = "variables6";
    variables[7] = "variables7";
    variables[8] = "variables8";
    variables[9] = "variables9";
    variables[10] = "variables10";
    variables[11] = "variables11";
    variables[12] = "variables12";
    variables[13] = "variables13";
    variables[14] = "variables14";
    variables[15] = "variables15";
    variables[16] = "variables16";
    variables[17] = "variables17";
    variables[18] = "variables18";
    variables[19] = "variables19";
    variables[20] = "variables20";
    variables[21] = "variables21";
    variables[22] = "variables22";
    averagedVariables[0] = "averagedVariables0";
    averagedVariables[1] = "averagedVariables1";
    averagedVariables[2] = "averagedVariables2";
    averagedVariables[3] = "averagedVariables3";
    averagedVariables[4] = "averagedVariables4";
    averagedVariables[5] = "averagedVariables5";
    averagedVariables[6] = "averagedVariables6";
    averagedVariables[7] = "averagedVariables7";
    averagedVariables[8] = "averagedVariables8";
    dt1Variables[0] = "dt1Variables0";
    dt1Variables[1] = "dt1Variables1";
    dt1Variables[2] = "dt1Variables2";
    dt1Variables[3] = "dt1Variables3";
    dt1Variables[4] = "dt1Variables4";
    dt1Variables[5] = "dt1Variables5";
    dt2Variables[0] = "dt2Variables0";
    dt2Variables[1] = "dt2Variables1";
    dt2Variables[2] = "dt2Variables2";
    dt2Variables[3] = "dt2Variables3";
    dt2Variables[4] = "dt2Variables4";
    dt2Variables[5] = "dt2Variables5";
    distributions[0] = "distributions0";
    distributions[1] = "distributions1";
    distributions[2] = "distributions2";
    distributions[3] = "distributions3";
    distributions[4] = "distributions4";
    distributions[5] = "distributions5";
    distributions[6] = "distributions6";
    distributions[7] = "distributions7";
    distributions[8] = "distributions8";
    distributions[9] = "distributions9";
    distributions[10] = "distributions10";
    distributions[11] = "distributions11";
    distributions[12] = "distributions12";
    distributions[13] = "distributions13";
    distributions[14] = "distributions14";
    distributions[15] = "distributions15";
    distributions[16] = "distributions16";
    distributions[17] = "distributions17";
    distributions[18] = "distributions18";
    distributions[19] = "distributions19";
    distributions[20] = "distributions20";
    distributions[21] = "distributions21";
    distributions[22] = "distributions22";
    distributions[23] = "distributions23";
    distributions[24] = "distributions24";
    distributions[25] = "distributions25";
    distributions[26] = "distributions26";
    oldDistributions[0] = "oldDistributions0";
    oldDistributions[1] = "oldDistributions1";
    oldDistributions[2] = "oldDistributions2";
    oldDistributions[3] = "oldDistributions3";
    oldDistributions[4] = "oldDistributions4";
    oldDistributions[5] = "oldDistributions5";
    oldDistributions[6] = "oldDistributions6";
    oldDistributions[7] = "oldDistributions7";
    oldDistributions[8] = "oldDistributions8";
    oldDistributions[9] = "oldDistributions9";
    oldDistributions[10] = "oldDistributions10";
    oldDistributions[11] = "oldDistributions11";
    oldDistributions[12] = "oldDistributions12";
    oldDistributions[13] = "oldDistributions13";
    oldDistributions[14] = "oldDistributions14";
    oldDistributions[15] = "oldDistributions15";
    oldDistributions[16] = "oldDistributions16";
    oldDistributions[17] = "oldDistributions17";
    oldDistributions[18] = "oldDistributions18";
    oldDistributions[19] = "oldDistributions19";
    oldDistributions[20] = "oldDistributions20";
    oldDistributions[21] = "oldDistributions21";
    oldDistributions[22] = "oldDistributions22";
    oldDistributions[23] = "oldDistributions23";
    oldDistributions[24] = "oldDistributions24";
    oldDistributions[25] = "oldDistributions25";
    oldDistributions[26] = "oldDistributions26";
    distributionsThermal[0] = "distributionsThermal0";
    distributionsThermal[1] = "distributionsThermal1";
    distributionsThermal[2] = "distributionsThermal2";
    distributionsThermal[3] = "distributionsThermal3";
    distributionsThermal[4] = "distributionsThermal4";
    distributionsThermal[5] = "distributionsThermal5";
    distributionsThermal[6] = "distributionsThermal6";
    distributionsThermal[7] = "distributionsThermal7";
    distributionsThermal[8] = "distributionsThermal8";
    distributionsThermal[9] = "distributionsThermal9";
    distributionsThermal[10] = "distributionsThermal10";
    distributionsThermal[11] = "distributionsThermal11";
    distributionsThermal[12] = "distributionsThermal12";
    distributionsThermal[13] = "distributionsThermal13";
    distributionsThermal[14] = "distributionsThermal14";
    distributionsThermal[15] = "distributionsThermal15";
    distributionsThermal[16] = "distributionsThermal16";
    distributionsThermal[17] = "distributionsThermal17";
    distributionsThermal[18] = "distributionsThermal18";
    distributionsThermal[19] = "distributionsThermal19";
    distributionsThermal[20] = "distributionsThermal20";
    distributionsThermal[21] = "distributionsThermal21";
    distributionsThermal[22] = "distributionsThermal22";
    distributionsThermal[23] = "distributionsThermal23";
    distributionsThermal[24] = "distributionsThermal24";
    distributionsThermal[25] = "distributionsThermal25";
    distributionsThermal[26] = "distributionsThermal26";
    oldDistributionsThermal[0] = "oldDistributionsThermal0";
    oldDistributionsThermal[1] = "oldDistributionsThermal1";
    oldDistributionsThermal[2] = "oldDistributionsThermal2";
    oldDistributionsThermal[3] = "oldDistributionsThermal3";
    oldDistributionsThermal[4] = "oldDistributionsThermal4";
    oldDistributionsThermal[5] = "oldDistributionsThermal5";
    oldDistributionsThermal[6] = "oldDistributionsThermal6";
    oldDistributionsThermal[7] = "oldDistributionsThermal7";
    oldDistributionsThermal[8] = "oldDistributionsThermal8";
    oldDistributionsThermal[9] = "oldDistributionsThermal9";
    oldDistributionsThermal[10] = "oldDistributionsThermal10";
    oldDistributionsThermal[11] = "oldDistributionsThermal11";
    oldDistributionsThermal[12] = "oldDistributionsThermal12";
    oldDistributionsThermal[13] = "oldDistributionsThermal13";
    oldDistributionsThermal[14] = "oldDistributionsThermal14";
    oldDistributionsThermal[15] = "oldDistributionsThermal15";
    oldDistributionsThermal[16] = "oldDistributionsThermal16";
    oldDistributionsThermal[17] = "oldDistributionsThermal17";
    oldDistributionsThermal[18] = "oldDistributionsThermal18";
    oldDistributionsThermal[19] = "oldDistributionsThermal19";
    oldDistributionsThermal[20] = "oldDistributionsThermal20";
    oldDistributionsThermal[21] = "oldDistributionsThermal21";
    oldDistributionsThermal[22] = "oldDistributionsThermal22";
    oldDistributionsThermal[23] = "oldDistributionsThermal23";
    oldDistributionsThermal[24] = "oldDistributionsThermal24";
    oldDistributionsThermal[25] = "oldDistributionsThermal25";
    oldDistributionsThermal[26] = "oldDistributionsThermal26";
  }

 public:
  const MChar* nDim;
  const MChar* noCells;
  const MChar* noConsistentCells;
  const MChar* domainId;
  const MChar* adjacentDomains;
  const MChar* windowCellCode;
  const MChar* originalDomain;
  const MChar* currentDomain;
  const MChar* originalId;
  const MChar* periodicDomain;
  const MChar* quasiPeriodicDomain;
  const MChar* noChildIds;
  const MChar* noCellsDirections;
  const MChar* noNghbrIds;
  const MChar* lengthLevel0;
  const MChar* maxLevel;
  const MChar* parentIds;
  const MChar* noCellChildIds;
  const MChar* childIds;
  const MChar* noCellNghbrIds;
  const MChar* nghbrIds;
  const MChar* intCoordinates;
  const MChar* level;
  const MChar* per;
  const MChar* qper;

  const MChar* noBndCnds;      // dimension
  const MChar* bndCndDim;      // dimension
  const MChar* bndCndIds;      // variable
  const MChar* bndCellIds;     // variable
  const MChar* bndCellSideIds; // variable

  const MChar* coordinates;
  const MChar* bodyCoordinates;

  const MChar* cellIds;
  const MChar* volume;
  const MChar* adaptation;
  const MChar* cutSideIds;
  const MChar* noCutSideIds;
  const MChar* nonFluidSideIds;
  const MChar* noNonFluidSideIds;
  const MChar* bndryCndId;
  const MChar* slope;
  const MChar* curvature;
  const MChar* area;
  const MChar* normalVctr;
  const MChar* cutCoordinates;
  const MChar* noCutPoints;
  const MChar* cutWedgeIds;
  const MChar* segmentIds;
  const MChar* bodyIds;
  const MChar* srfcs;
  const MChar* noSrfcs;

  const MChar* rightHandSide;
  const MChar* rcnstrctnConstants;
  const MChar* noRcnstrctnNghbrIds;
  const MChar* rcnstrctnNghbrIds;
  const MChar* u_Velocity;
  const MChar* v_Velocity;
  const MChar* w_Velocity;
  const MChar* temperature;
  const MChar* energy;
  const MChar* rho;
  const MChar* passiveScalar;
  const MChar* minValues;
  const MChar* maxValues;
  const MChar* noSamples;
  const MChar* globalTimeStep;
  const MChar* del;
  const MChar* time;
  const MChar* physicalTime;
  const MChar* firstMaxResidual;
  const MChar* firstAvrgResidual;
  const MChar* levelSetFunction;
  const MChar* oldLevelSetFunction;
  const MChar* extensionVelocity;
  const MChar* semiLagrangeXShiftRef;
  const MChar* cellIsInactive;

  const MChar* flowCellId;
  const MChar* bndryId;

  const MChar* variables[23];
  const MChar* averagedVariables[9];
  const MChar* oldVariables[6];
  const MChar* dt1Variables[6];
  const MChar* dt2Variables[6];
  const MChar* distributions[27];
  const MChar* oldDistributions[27];
  const MChar* distributionsThermal[27];
  const MChar* oldDistributionsThermal[27];
};

#endif // CARTESIANNETCDF_H
